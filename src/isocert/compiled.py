"""Compiled envelope tables: compile a record once, answer every later query by table lookup.

U is convex in ln r, so linear interpolation of (ln r_j, U_j) over-estimates it; L satisfies
L(r) >= L(r_j) + q_lo(r_j) ln(r / r_j) for r >= r_j because every member is nondecreasing.  Both interpolants are
therefore outer bounds, and intervals computed from tables contain the exact certified intervals.

Outside the tabulated range the table falls back to bounds that hold for every member of the class:
  r < r_1:  U(r) <= min(H r, U(r_1)),            L(r) >= L(r_1) r / r_1        (q(t)/t is nonincreasing and <= H)
  r > r_m:  U(r) <= U(r_m) + q_hi(r_m)(r/r_m-1), L(r) >= L(r_m) + q_lo(r_m) ln(r/r_m)
so a CompiledRecord is self-contained: it can be saved as JSON and queried without the solver.
"""
import json
import math
import numpy as np

_SIG = 7


def _up(v):      # outward rounding to _SIG significant digits
    if v == 0 or not math.isfinite(v): return v
    s = 10 ** (_SIG - 1 - math.floor(math.log10(abs(v))))
    return math.ceil(v * s) / s


def _down(v):
    if v == 0 or not math.isfinite(v): return v
    s = 10 ** (_SIG - 1 - math.floor(math.log10(abs(v))))
    return math.floor(v * s) / s


class CompiledRecord:
    """Sound lookup tables for the envelopes U(r), L(r) of one record at one Henry cap."""

    def __init__(self, r, U, L, q_lo, H, q_hi_top, meta=None):
        self.r = np.asarray(r, float); self.lnr = np.log(self.r)
        self.Uj = np.asarray(U, float); self.Lj = np.asarray(L, float); self.qlo = np.asarray(q_lo, float)
        self.H = float(H); self.q_hi_top = float(q_hi_top); self.meta = dict(meta or {})

    @classmethod
    def from_source(cls, src, r_lo=None, r_hi=None, per_decade=12, meta=None):
        """Tabulate a solver object (isocert.envelopes.StSource).  Default range: 0.1*min(p_1, 0.005) to 1e4*P_max."""
        r_lo = min(src.p[0], 0.005) * 0.1 if r_lo is None else r_lo
        r_hi = src.P * src.RCAP if r_hi is None else r_hi
        n = max(2, int(math.ceil(math.log10(r_hi / r_lo) * per_decade)) + 1)
        r = np.array([float('%.9g' % x) for x in np.geomspace(r_lo, r_hi, n)])     # the stored nodes are the evaluated nodes
        return cls(r, [src.U(x) for x in r], [src.L(x) for x in r], [src.q_lo(x) for x in r], src.H, src.q_hi(r[-1]), meta)

    def U(self, r):
        if r <= self.r[0]:
            return min(self.H * r, self.Uj[0])
        if r >= self.r[-1]:
            return self.Uj[-1] + self.q_hi_top * (r / self.r[-1] - 1.0)
        j = int(np.searchsorted(self.r, r, side='right')) - 1
        t = (math.log(r) - self.lnr[j]) / (self.lnr[j + 1] - self.lnr[j])
        return (1 - t) * self.Uj[j] + t * self.Uj[j + 1]

    def L(self, r):
        if r <= self.r[0]:
            return self.Lj[0] * r / self.r[0]
        j = min(int(np.searchsorted(self.r, r, side='right')) - 1, len(self.r) - 1)
        return self.Lj[j] + self.qlo[j] * (math.log(r) - self.lnr[j])

    # ---- serialisation (values rounded outward, so a saved table is still a sound bound)
    def to_dict(self):
        return dict(r=[float('%.9g' % x) for x in self.r], U=[_up(float(x)) for x in self.Uj], L=[_down(float(x)) for x in self.Lj],
                    q_lo=[_down(float(x)) for x in self.qlo], H=_up(self.H), q_hi_top=_up(self.q_hi_top), meta=self.meta)

    @classmethod
    def from_dict(cls, d):
        return cls(d['r'], d['U'], d['L'], d['q_lo'], d['H'], d['q_hi_top'], d.get('meta'))

    def save(self, path):
        with open(path, 'w') as f: json.dump(self.to_dict(), f)

    @classmethod
    def load(cls, path):
        with open(path) as f: return cls.from_dict(json.load(f))
