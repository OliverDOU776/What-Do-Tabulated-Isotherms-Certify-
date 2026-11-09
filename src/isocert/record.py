"""High-level interface: a tabulated record, its class of consistent isotherms, and queries on pairs of records."""
import math
import numpy as np

from .envelopes import StSource, misfit_floor, min_henry
from .compiled import CompiledRecord
from .iast import interval, decide, lnS_bounds

EPS_FLOOR = 0.01      # tolerance floor: 1% of full scale
EPS_MULT = 3.0        # tolerance = max(EPS_MULT * misfit floor, EPS_FLOOR); calibrated on held-out data in the paper
EPS_ADMIT = 0.03      # a record is admitted when its misfit floor is at most 3% of full scale
MAX_POINTS = 40       # longer records are thinned to 40 points by index


def thin(p, q, nmax=MAX_POINTS):
    p = np.asarray(p, float); q = np.asarray(q, float)
    if len(p) <= nmax: return p, q
    idx = np.unique(np.round(np.linspace(0, len(p) - 1, nmax)).astype(int))
    return p[idx], q[idx]


class Record:
    """A tabulated pure-gas isotherm (pressures p, loadings q) and the class of isotherms consistent with it.

    The class contains every site mixture q(p) = int phi_b(p) dnu(b) with |q(p_i) - q_i| <= e_i and Henry constant
    at most H, where e_i = eps * max q, eps = max(multiplier * misfit floor, 1%), and H = kappa * K_min
    (K_min is the smallest Henry constant consistent with the table).

    Parameters
    ----------
    p, q : increasing positive pressures and nonnegative loadings (any consistent units).
    kernel : 'langmuir' (default; multi-site Langmuir, UNILAN, Toth, ...) or 'hinge' (shape only: concave, nondecreasing).
    kappa : Henry-cap multiplier (default 10).  Pass H to fix the cap in absolute units instead.
    eps : full-scale tolerance; by default derived from the misfit floor.
    """

    def __init__(self, p, q, *, kernel='langmuir', kappa=10.0, H=None, eps=None, multiplier=EPS_MULT, name=None):
        p = np.asarray(p, float); q = np.asarray(q, float)
        if p.ndim != 1 or len(p) != len(q) or len(p) < 3: raise ValueError('need at least three (p, q) points')
        if np.any(np.diff(p) <= 0) or p[0] <= 0: raise ValueError('pressures must be positive and strictly increasing')
        if np.any(q < 0) or q.max() <= 0: raise ValueError('loadings must be nonnegative and not all zero')
        self.p_all, self.q_all = p, q
        self.p, self.q = thin(p, q)
        self.kernel, self.kappa, self.name = kernel, float(kappa), name
        self.eps_star = misfit_floor(self.p, self.q, kernel=kernel)
        self.admitted = bool(math.isfinite(self.eps_star) and self.eps_star <= EPS_ADMIT)
        self.eps = float(eps) if eps is not None else max(multiplier * self.eps_star, EPS_FLOOR)
        self._H = H
        self._env = None

    @property
    def envelope(self):
        """Solver object exposing U(r), L(r), q_lo(r), q_hi(r) (built on first use)."""
        if self._env is None:
            self._env = (StSource(self.p, self.q, self.eps, H=self._H, kernel=self.kernel) if self._H is not None
                         else StSource(self.p, self.q, self.eps, kappa=self.kappa, kernel=self.kernel))
            if not self._env.feasible: raise ValueError('no isotherm of the class fits the record within tolerance and cap')
        return self._env

    @property
    def H(self): return self.envelope.H

    @property
    def K_min(self):
        return self.envelope.Kmin if self.envelope.Kmin is not None else min_henry(self.p, self.q, self.eps, kernel=self.kernel)

    def U(self, r): return self.envelope.U(r)
    def L(self, r): return self.envelope.L(r)

    def band(self, pressures):
        """Certified loading band: arrays (q_lo, q_hi) such that every isotherm of the class lies between them."""
        pr = np.atleast_1d(np.asarray(pressures, float)); e = self.envelope
        return np.array([e.q_lo(x) for x in pr]), np.array([e.q_hi(x) for x in pr])

    def sites(self, r, which='U'):
        """Sites (affinity, capacity, Henry mass) of an extremal isotherm attaining the envelope at r."""
        return self.envelope.extremal_sites(r, which)

    def compile(self, per_decade=12):
        """Tabulate the envelopes once; the result answers later queries by lookup and can be saved as JSON."""
        return CompiledRecord.from_source(self.envelope, per_decade=per_decade,
                                          meta=dict(name=self.name, kernel=self.kernel, kappa=self.kappa, eps=self.eps,
                                                    eps_star=self.eps_star, p_min=float(self.p[0]), p_max=float(self.p[-1])))

    def with_cap(self, kappa):
        """The same record under another Henry-cap multiplier."""
        return Record(self.p_all, self.q_all, kernel=self.kernel, kappa=kappa, eps=self.eps, name=self.name)

    def with_point(self, p_new, q_new):
        """The record after one more measurement.  Absolute tolerance and Henry cap are kept, so the class can only shrink."""
        i = int(np.searchsorted(self.p, p_new))
        p2 = np.insert(self.p, i, p_new); q2 = np.insert(self.q, i, q_new)
        r = Record.__new__(Record); r.__dict__.update(self.__dict__)
        r.p_all, r.q_all, r.p, r.q = p2, q2, p2, q2
        r._H = self.H
        r._env = StSource(p2, q2, self.eps, H=self.H, e=np.full(len(p2), self.eps * self.q.max()), kernel=self.kernel)
        return r


def critical_cap(A, B, P, y, tau, kappa_max=1e5, steps=8):
    """Largest Henry-cap multiplier (up to kappa_max) at which the answer certified at the records' own caps is still
    certified.  Returns None when the query is UNDETERMINED at the current caps.  Two envelope LPs per cap tested."""
    base = decide(A, B, P, y, tau)
    if base == 'UNDETERMINED': return None
    holds = lambda k: decide(A.with_cap(k), B.with_cap(k), P, y, tau) == base
    k0 = max(A.kappa, B.kappa)
    if holds(kappa_max): return float(kappa_max)
    lo, hi = math.log(k0), math.log(kappa_max)
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if holds(math.exp(mid)): lo = mid
        else: hi = mid
    return math.exp(lo)


def next_measurement(A, B, P, y, candidates=None, top=None):
    """One-step lookahead: which single new measurement (gas, pressure) would narrow the certified interval most?

    For every candidate pressure of each gas, the new loading is posited at the middle of the current certified band
    and the interval is recomputed.  Returns dicts sorted by predicted width of ln S:
    {'gas': 'A' | 'B', 'pressure', 'posited_loading', 'S_lo', 'S_hi', 'width'}.  A myopic heuristic, not an optimal plan."""
    cands = np.geomspace(1e-4, 100.0, 25) if candidates is None else np.asarray(candidates, float)
    out = []
    for g, (R, other) in (('A', (A, B)), ('B', (B, A))):
        for pc in cands:
            if np.any(np.abs(np.log(R.p / pc)) < 0.05): continue        # already measured there
            lo, hi = R.band([pc]); qm = 0.5 * float(lo[0] + hi[0])
            try:
                R2 = R.with_point(pc, qm)
                iv = interval(R2, other, P, y) if g == 'A' else interval(other, R2, P, y)
            except Exception:
                continue
            out.append(dict(gas=g, pressure=float(pc), posited_loading=qm, S_lo=iv.S_lo, S_hi=iv.S_hi, width=iv.width))
    out.sort(key=lambda d: d['width'])
    return out[:top] if top else out
