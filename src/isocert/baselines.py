"""Baselines that current practice (or a careful version of it) would report.

* Six parametric families fitted by least squares: Langmuir, dual-site Langmuir, Sips, Toth, UNILAN, Freundlich.
  `fit_all` returns every fit with its full-scale RMSE and maximum residual; an "ensemble" is the set of fits
  whose maximum residual is within the record's tolerance.
* `DualSiteSet`: a grid enumeration of ALL dual-site Langmuir isotherms within the same tolerance and Henry cap
  (an inner approximation without a certificate; the strongest same-information baseline in the paper).
"""
import math
import numpy as np
from scipy.optimize import least_squares

_GLX, _GLW = np.polynomial.legendre.leggauss(96)


def _pi_numeric(qfun, r, a0=1.0):
    """int_0^r q(t)/t dt.  Prefix [0, t0] handled with q ~ C t^a0 (a0 = small-p exponent);
    the rest by composite Gauss-Legendre in u = ln t."""
    t0 = r * 1e-12
    pref = qfun(np.array([t0]))[0] / a0
    lo, hi = math.log(t0), math.log(r)
    tot = 0.0
    nseg = 6
    edges = np.linspace(lo, hi, nseg + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        u = 0.5 * (b - a) * _GLX + 0.5 * (b + a)
        tot += 0.5 * (b - a) * np.dot(_GLW, qfun(np.exp(u)))
    return pref + tot


class Model:
    name = ''; npar = 0; in_class = True
    def __init__(self, theta): self.theta = np.asarray(theta, float)
    def henry(self): return float('nan')


class Langmuir(Model):
    name = 'Langmuir'; npar = 2
    def q(self, p): qs, b = np.exp(self.theta); return qs * b * p / (1 + b * p)
    def Pi(self, r): qs, b = np.exp(self.theta); return qs * math.log1p(b * r)
    def henry(self): qs, b = np.exp(self.theta); return qs * b
    @staticmethod
    def starts(p, q): return [[math.log(q.max() * f), math.log(g / p[-1])] for f in (1.2, 3.0) for g in (0.3, 3.0, 30.0)]


class DSL(Model):
    name = 'DSL'; npar = 4
    def q(self, p): q1, b1, q2, b2 = np.exp(self.theta); return q1 * b1 * p / (1 + b1 * p) + q2 * b2 * p / (1 + b2 * p)
    def Pi(self, r): q1, b1, q2, b2 = np.exp(self.theta); return q1 * math.log1p(b1 * r) + q2 * math.log1p(b2 * r)
    def henry(self): q1, b1, q2, b2 = np.exp(self.theta); return q1 * b1 + q2 * b2
    @staticmethod
    def starts(p, q):
        qm = q.max(); P = p[-1]
        return [[math.log(qm * a), math.log(g1 / P), math.log(qm * c), math.log(g2 / P)]
                for (a, c) in ((0.5, 1.5), (1.0, 3.0)) for (g1, g2) in ((30.0, 0.3), (300.0, 1.0), (5.0, 0.05))]


class Sips(Model):
    name = 'Sips'; npar = 3; in_class = True
    def q(self, p): qs, b, n = np.exp(self.theta[0]), np.exp(self.theta[1]), 1 / (1 + np.exp(-self.theta[2])); x = (b * p) ** n; return qs * x / (1 + x)
    def Pi(self, r): qs, b, n = np.exp(self.theta[0]), np.exp(self.theta[1]), 1 / (1 + np.exp(-self.theta[2])); return qs / n * math.log1p((b * r) ** n)
    def henry(self): n = 1 / (1 + np.exp(-self.theta[2])); return float('inf') if n < 1.0 else np.exp(self.theta[0] + self.theta[1])   # finite only when the exponent is exactly 1 (Langmuir)
    @staticmethod
    def starts(p, q): return [[math.log(q.max() * f), math.log(g / p[-1]), s] for f in (1.5, 4.0) for g in (0.3, 3.0) for s in (0.5, 3.0)]


class Toth(Model):
    name = 'Toth'; npar = 3
    def _par(self): return np.exp(self.theta[0]), np.exp(self.theta[1]), 1 / (1 + np.exp(-self.theta[2]))
    def q(self, p): qs, b, t = self._par(); return qs * b * p / (1 + (b * p) ** t) ** (1 / t)
    def Pi(self, r): return _pi_numeric(self.q, r, 1.0)
    def henry(self): qs, b, t = self._par(); return qs * b
    @staticmethod
    def starts(p, q): return [[math.log(q.max() * f), math.log(g / p[-1]), s] for f in (1.5, 4.0) for g in (1.0, 30.0) for s in (0.0, 2.0)]


class Unilan(Model):
    name = 'UNILAN'; npar = 3
    def _par(self): return np.exp(self.theta[0]), np.exp(self.theta[1]), np.exp(self.theta[2])
    def q(self, p):
        qs, b, s = self._par(); s = min(max(s, 1e-3), 30.0)
        return qs / (2 * s) * (np.log1p(b * math.exp(s) * p) - np.log1p(b * math.exp(-s) * p))
    def Pi(self, r): return _pi_numeric(self.q, r, 1.0)
    def henry(self): qs, b, s = self._par(); s = min(max(s, 1e-3), 30.0); return qs * b * math.sinh(s) / s
    @staticmethod
    def starts(p, q): return [[math.log(q.max() * f), math.log(g / p[-1]), math.log(s)] for f in (1.5, 4.0) for g in (0.3, 3.0) for s in (1.0, 4.0)]


class Freundlich(Model):
    name = 'Freundlich'; npar = 2
    def _par(self): return np.exp(self.theta[0]), 1 / (1 + np.exp(-self.theta[1]))
    def q(self, p): k, a = self._par(); return k * p ** a
    def Pi(self, r): k, a = self._par(); return k * r ** a / a
    def henry(self): return float('inf')
    @staticmethod
    def starts(p, q): return [[math.log(q[-1] / p[-1] ** a), math.log(a / (1 - a))] for a in (0.3, 0.6, 0.9)]


FAMILIES = [Langmuir, DSL, Sips, Toth, Unilan, Freundlich]


def fit(cls, p, q):
    """Least-squares fit (absolute residuals); returns (model, rmse_fs, sup_fs)."""
    p = np.asarray(p, float); q = np.asarray(q, float); qm = q.max()
    best = None
    for x0 in cls.starts(p, q):
        try:
            res = least_squares(lambda th: (cls(th).q(p) - q) / qm, np.array(x0, float), method='trf',
                                max_nfev=300, x_scale=1.0, bounds=(-60, 60))
        except Exception:
            continue
        if best is None or res.cost < best.cost:
            best = res
    if best is None:
        return None, float('nan'), float('nan')
    m = cls(best.x)
    r = (m.q(p) - q) / qm
    return m, float(np.sqrt(np.mean(r ** 2))), float(np.max(np.abs(r)))


def fit_all(p, q):
    out = {}
    for cls in FAMILIES:
        m, rmse, sup = fit(cls, p, q)
        if m is not None and np.all(np.isfinite(m.theta)):
            out[cls.name] = (m, rmse, sup)
    return out


# ---------------------------------------------------------------------------------------------------
# Dual-site Langmuir parameter set: q(p) = m1 p/(1+b1 p) + m2 p/(1+b2 p), m1 + m2 <= H, |q(p_i) - q_i| <= e_i.
# Variables are scaled as in the envelope solver (u_j = loading contributed by site j at P_max, over max q).
def _a(x, P, b):          # phi_b(x)/phi_b(P) for pressures x (rows) and affinities b (cols)
    x = np.atleast_1d(np.asarray(x, float)); return 1.0 - (1.0 - x / P)[:, None] / (1.0 + np.outer(x, b))


def _c(r, P, b):          # psi_r(b)/phi_b(P)
    v = b * r; ell = np.where(v > 1e-8, np.log1p(v) / np.maximum(v, 1e-300), 1.0 - 0.5 * v)
    return (r / P) * ell + np.log1p(v)


class DSLSet:
    def __init__(self, p, q, eps, H, per_decade=10, nu=40, nb2=44, nu2=56):
        self.p = np.asarray(p, float); self.q = np.asarray(q, float); self.P = float(self.p[-1]); self.qm = float(self.q.max())
        e = eps * self.qm; self.up = (self.q + e) / self.qm; self.lo = (self.q - e) / self.qm; self.Hn = H / self.qm
        lo_d = math.log10(1e-4 / self.P); hi_d = math.log10(1e4 / self.p[0])
        g = np.concatenate([[0.0], np.logspace(lo_d, hi_d, int(round((hi_d - lo_d) * per_decade)) + 1)])
        parts = [self._search(g, g, nu, tri=True)]
        if len(parts[0][0]):                      # second level: refine the bounding box of feasible affinity pairs
            B1, B2 = parts[0][0], parts[0][1]; st = 10 ** (1.0 / per_decade)
            def rng(B, n):
                pos = B[B > 0]
                if len(pos) == 0: return np.array([0.0])
                gg = np.logspace(math.log10(pos.min() / st), math.log10(pos.max() * st), n)
                return np.concatenate([[0.0], gg]) if (B == 0).any() or pos.min() <= g[1] * 1.0001 else gg
            parts.append(self._search(rng(B1, nb2), rng(B2, nb2), nu2, tri=False))
        self.B1, self.B2, self.U1, self.U2lo, self.U2hi = [np.concatenate([x[k] for x in parts]) for k in range(5)]
        self.n = len(self.B1); self.empty = self.n == 0

    def _search(self, g1, g2, nu, tri):
        P, Hn = self.P, self.Hn
        A1 = _a(self.p, P, g1); A2 = _a(self.p, P, g2); M1 = (1 + g1 * P) / (Hn * P); M2 = (1 + g2 * P) / (Hn * P)
        out = [[], [], [], [], []]
        for j in range(len(g1)):
            a1 = A1[:, j]
            u1max = min(float((self.up / a1).min()), 1.0 / M1[j])
            if u1max <= 0: continue
            u1 = np.linspace(0.0, u1max, nu)
            ru = self.up[:, None] - a1[:, None] * u1[None, :]; rl = self.lo[:, None] - a1[:, None] * u1[None, :]      # n x K
            hi = (ru[:, :, None] / A2[:, None, :]).min(0); hi = np.minimum(hi, (1.0 - u1 * M1[j])[:, None] / M2[None, :])
            lo = np.maximum((rl[:, :, None] / A2[:, None, :]).max(0), 0.0)
            feas = lo <= hi
            if tri: feas[:, :j] = False                # unordered pairs once
            k, j2 = np.where(feas)
            if len(k) == 0: continue
            out[0].append(np.full(len(k), g1[j])); out[1].append(g2[j2]); out[2].append(u1[k]); out[3].append(lo[k, j2]); out[4].append(hi[k, j2])
        return [np.concatenate(x) if x else np.array([]) for x in out]

    def band(self, r):
        """(q_lo, q_hi) of the parameter set at pressures r."""
        r = np.atleast_1d(np.asarray(r, float)); lo = np.empty(len(r)); hi = np.empty(len(r))
        for i, x in enumerate(r):
            a1 = _a([x], self.P, self.B1)[0]; a2 = _a([x], self.P, self.B2)[0]
            lo[i] = (self.U1 * a1 + self.U2lo * a2).min(); hi[i] = (self.U1 * a1 + self.U2hi * a2).max()
        return lo * self.qm, hi * self.qm

    def U(self, r):
        return float((self.U1 * _c(r, self.P, self.B1) + self.U2hi * _c(r, self.P, self.B2)).max()) * self.qm

    def L(self, r):
        return float((self.U1 * _c(r, self.P, self.B1) + self.U2lo * _c(r, self.P, self.B2)).min()) * self.qm


DualSiteSet = DSLSet
