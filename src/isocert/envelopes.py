"""Certified spreading-pressure and loading envelopes of site-mixture classes.

q(p)  = int_{[0,inf)} phi_b(p) dnu(b),   nu >= 0,   nu([0,inf)) = Henry constant <= H
Pi(r) = int_{[0,inf)} psi_r(b) dnu(b)                (reduced spreading pressure)
Data:   |q(p_i) - q_i| <= e_i                         (full-scale tolerance e_i = eps * max_i q_i)

Kernels: 'langmuir' (phi_b = p/(1+bp): Stieltjes class, contains multi-site Langmuir, UNILAN, Toth) and
'hinge' (phi_b = min(p, 1/b): all concave nondecreasing isotherms through the origin).

Each envelope value is a generalised moment problem.  It is solved as a linear programme on an affinity grid with
column generation, and the returned number is the dual objective minus a certificate margin that holds for EVERY
measure on [0, inf): the reduced cost is bounded below on the whole half-line by interval bounds with adaptive
bisection (class _RC), including both tails.  Upper (lower) envelope values are therefore upper (lower) bounds for
the continuous class, up to double-precision rounding (an outward relative slack of 2e-14 is included).
The numerical kernel is the one used for the paper's experiments.
"""
import math
import numpy as np
from scipy.optimize import linprog


def _grid(p1, Pmax, kappa, per_decade=16, lo_dec=5.0):
    lo = math.log10(1.0 / Pmax) - lo_dec
    hi = math.log10(1.0 / p1) + 3.0 + max(math.log10(max(kappa, 1.0)), 0.0)
    n = int(round((hi - lo) * per_decade)) + 1
    return np.concatenate([[0.0], np.logspace(lo, hi, n)])


def _phi(p, b, kernel='langmuir'):
    """Unit-Henry-mass site response at pressures p (rows) for affinities b (cols).
    langmuir: p/(1+bp);  hinge: min(p, 1/b)  (concave class = mixtures of hinge sites)."""
    if kernel == 'langmuir':
        return p[:, None] / (1.0 + np.outer(p, b))
    with np.errstate(divide='ignore'):
        s = np.where(b > 0, 1.0 / np.maximum(b, 1e-300), np.inf)
    return np.minimum(p[:, None], s[None, :])


def _kern(r, b, kernel='langmuir'):
    """Spreading pressure of a unit-Henry-mass site at horizon r."""
    out = np.full(b.shape, float(r))
    m = b > 0
    if kernel == 'langmuir':
        out[m] = np.log1p(b[m] * r) / b[m]
    else:
        br = b[m] * r
        out[m] = np.where(br <= 1.0, r, (1.0 + np.log(np.maximum(br, 1.0))) / b[m])
    return out


def _qk(r, b, kernel='langmuir'):
    return _phi(np.array([float(r)]), b, kernel)[0]


def misfit_floor(p, q, per_decade=16, kernel='langmuir'):
    """eps* = smallest full-scale tolerance for which SOME site distribution matches the data."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    b = _grid(p[0], p[-1], 1e4, per_decade)
    if kernel == 'hinge': b = np.unique(np.concatenate([b, 1.0 / p]))
    qm = q.max(); s = qm / p[-1]
    A = _phi(p, b, kernel) * np.maximum(b, 1.0 / p[-1])[None, :] if kernel == 'langmuir' else _phi(p, b, kernel) / p[-1]
    A[:, 0] = p / p[-1]
    m = A.shape[1]
    c = np.zeros(m + 1); c[-1] = 1.0
    one = np.ones((len(p), 1))
    A_ub = np.vstack([np.hstack([A, -one]), np.hstack([-A, -one])])
    b_ub = np.concatenate([q / qm, -q / qm])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method='highs')
    return float(res.x[-1]) if res.status == 0 else float('nan')


def min_henry(p, q, eps, per_decade=16, kernel='langmuir'):
    """Smallest Henry constant consistent with the data at full-scale tolerance eps."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    b = _grid(p[0], p[-1], 1e4, per_decade)
    if kernel == 'hinge': b = np.unique(np.concatenate([b, 1.0 / p]))
    qm = q.max(); e = eps * qm
    s = qm / p[-1]                                   # Henry scale
    A = _phi(p, b, kernel) * s / qm
    A_ub = np.vstack([A, -A]); b_ub = np.concatenate([(q + e) / qm, -(q - e) / qm])
    res = linprog(np.ones(len(b)), A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method='highs')
    if res.status != 0:
        return float('nan')
    return float(res.fun) * s


def _ell(v):
    """ln(1+v)/v, with the limits 1 at v = 0 and 0 at v = inf (decreasing in v)."""
    v = np.asarray(v, float); out = np.ones_like(v)
    fin = np.isfinite(v); big = fin & (v > 1e-8)
    out[big] = np.log1p(v[big]) / v[big]
    sm = fin & ~big; out[sm] = 1.0 - 0.5 * v[sm]
    out[~fin] = 0.0
    return out


def _f4(u):
    """u/(1+u)^2: increasing on [0,1], decreasing on [1,inf), maximum 1/4."""
    u = np.asarray(u, float); out = np.zeros_like(u); fin = np.isfinite(u)
    out[fin] = u[fin] / (1.0 + u[fin]) ** 2
    return out


def _f4_range(ulo, uhi):
    a = _f4(ulo); b = _f4(uhi)
    return np.minimum(a, b), np.where((ulo <= 1.0) & (uhi >= 1.0), 0.25, np.maximum(a, b))


class _RC:
    """Reduced cost of the site with affinity b, in the scaled variables of StSource._cols,
        rc(b) = s*c(b) + sum_i w_i a_i(b) + lm*m(b),   a_i = phi_b(p_i)/phi_b(P),  m = 1/(Hn phi_b(P)),  c = g_r(b)/phi_b(P),
    with rigorous lower bounds on intervals of b.  Facts used (both kernels): a_i and m are nondecreasing in b;
    Langmuir: a_i = 1-(1-rho_i)/(1+b p_i), m = (1+bP)/(Hn P), psi-objective c = (r/P) ell(br) + ln(1+br) (decreasing + increasing part),
              phi-objective c = 1-(1-r/P)/(1+br) (monotone);
    hinge:    a_i = rho_i (b<=1/P), p_i b (1/P<=b<=1/p_i), 1 (b>=1/p_i); m = max(1/P,b)/Hn;
              psi-objective c = psi_r(b) max(1/P,b) with psi_r nonincreasing; phi-objective c = min(r,1/b) max(1/P,b)."""

    def __init__(self, src, kind, r, s, lam):
        self.S = src; self.kind = kind; self.r = float(r); self.s = float(s); self.kernel = src.kernel
        n = src.n; w = lam[:n] - lam[n:2 * n]
        self.w = w; self.wp = np.maximum(w, 0.0); self.wn = np.minimum(w, 0.0); self.lm = float(lam[2 * n])
        self.p = src.p; self.P = float(src.P); self.Hn = float(src.Hn); self.rho = self.p / self.P
        self.wabs = float(np.abs(w).sum())

    # ---- elementary functions (b: 1-D array, 0 <= b <= inf)
    def A(self, b):
        if self.kernel == 'langmuir':
            return 1.0 - (1.0 - self.rho)[:, None] / (1.0 + np.outer(self.p, b))
        return np.where(b[None, :] <= 1.0 / self.P, self.rho[:, None], np.minimum(np.outer(self.p, b), 1.0))

    def m(self, b):
        if self.kernel == 'langmuir': return (1.0 + b * self.P) / (self.Hn * self.P)
        return np.maximum(b, 1.0 / self.P) / self.Hn

    def _psi_h(self, b):                 # hinge psi_r(b): r if br <= 1 else (1+ln(br))/b ; 0 at b = inf
        out = np.full(b.shape, self.r); br = b * self.r
        mk = np.isfinite(b) & (br > 1.0); out[mk] = (1.0 + np.log(br[mk])) / b[mk]; out[~np.isfinite(b)] = 0.0
        return out

    def c(self, b):
        r, P = self.r, self.P
        if self.kernel == 'langmuir':
            v = b * r
            return (r / P) * _ell(v) + np.log1p(v) if self.kind == 'psi' else 1.0 - (1.0 - r / P) / (1.0 + v)
        if self.kind == 'psi':
            br = b * r
            with np.errstate(divide='ignore', invalid='ignore'):
                hi = (1.0 + np.log(np.maximum(br, 1.0))) * np.maximum(1.0 / np.maximum(b * P, 1e-300), 1.0)
            return np.where(br <= 1.0, r * np.maximum(1.0 / P, b), hi)
        with np.errstate(divide='ignore'):
            inv = np.where(b > 0, 1.0 / np.maximum(b, 1e-300), np.inf)
        return np.where(b <= 1.0 / P, np.minimum(r, inv) / P, np.minimum(r * b, 1.0))

    def c_bounds(self, blo, bhi):
        """(min, max) of c on [blo, bhi] from monotone pieces."""
        r, P = self.r, self.P
        if self.kernel == 'langmuir':
            vlo, vhi = blo * r, bhi * r
            if self.kind == 'psi':
                return (r / P) * _ell(vhi) + np.log1p(vlo), (r / P) * _ell(vlo) + np.log1p(vhi)
            a, b_ = self.c(blo), self.c(bhi)
            return np.minimum(a, b_), np.maximum(a, b_)
        k_lo, k_hi = np.maximum(1.0 / P, blo), np.maximum(1.0 / P, bhi)
        if self.kind == 'psi':
            return self._psi_h(bhi) * k_lo, self._psi_h(blo) * k_hi
        with np.errstate(divide='ignore'):
            ilo = np.where(blo > 0, 1.0 / np.maximum(blo, 1e-300), np.inf); ihi = np.where(bhi > 0, 1.0 / np.maximum(bhi, 1e-300), np.inf)
        return np.minimum(r, ihi) * k_lo, np.minimum(np.minimum(r, ilo) * k_hi, max(r / P, 1.0))   # c <= max(r/P, 1) always

    def rc(self, b):
        return self.s * self.c(b) + self.w @ self.A(b) + self.lm * self.m(b)

    def _slack(self, cabs, mhi):
        return 2e-14 * (cabs + self.wabs + self.lm * mhi + 1.0)

    # ---- rigorous bounds
    def mono(self, blo, bhi):
        """Lower bound of rc on [blo, bhi] (bhi finite) from monotonicity of every term."""
        cmin, cmax = self.c_bounds(blo, bhi)
        sc = cmin if self.s > 0 else -cmax
        return sc + self.wp @ self.A(blo) + self.wn @ self.A(bhi) + self.lm * self.m(blo) - self._slack(np.maximum(np.abs(cmin), np.abs(cmax)), self.m(bhi))

    def tail(self, B):
        """Lower bound of rc on [B, inf)."""
        r, P, lm, Hn = self.r, self.P, self.lm, self.Hn
        Bv = np.array([float(B)])
        if self.kernel == 'langmuir':
            base = float(self.wp @ self.A(Bv)[:, 0] + self.wn.sum())
            if self.kind == 'psi' and self.s > 0:
                val = math.log1p(B * r) + base + lm * float(self.m(Bv)[0])
            elif self.kind == 'psi':
                if lm <= 0: return -math.inf
                bs = max(B, Hn / lm - 1.0 / r)          # minimiser of lm b/Hn - ln(1+br) on [B, inf)
                val = -(r / P) * float(_ell(np.array([B * r]))[0]) + base + lm / (Hn * P) + lm * bs / Hn - math.log1p(bs * r)
            else:
                cB = float(self.c(Bv)[0]); val = (min(cB, 1.0) if self.s > 0 else -max(cB, 1.0)) + base + lm * float(self.m(Bv)[0])
            return val - float(self._slack(abs(val), float(self.m(Bv)[0])))
        assert B >= max(1.0 / P, 1.0 / r, 1.0 / self.p[0])   # beyond every kink: a_i = 1, m = b/Hn, c = 1+ln(br) or 1
        sw = float(self.w.sum())
        if self.kind == 'psi' and self.s > 0:
            val = 1.0 + math.log(B * r) + sw + lm * B / Hn
        elif self.kind == 'psi':
            if lm <= 0: return -math.inf
            bs = max(B, Hn / lm); val = -1.0 - math.log(bs * r) + sw + lm * bs / Hn
        else:
            val = self.s + sw + lm * B / Hn
        return val - float(self._slack(abs(val), B / Hn))

    def centered(self, blo, bhi):
        """Lower bound of rc on each [blo, bhi] (0 < blo < bhi < inf): max of the monotone bound and the mean-value bound
        rc(bc) - (ln(bhi/blo)/2) * sup|d rc/d ln b|, the derivative being enclosed by interval arithmetic.
        For the hinge kernel every interval must lie inside one smooth piece (the node set contains all kinks)."""
        r, P, lm, Hn, p = self.r, self.P, self.lm, self.Hn, self.p
        bc = np.sqrt(blo * bhi); half = 0.5 * np.log(bhi / blo)
        Ac = self.A(bc); cc = self.c(bc); mc = self.m(bc)
        rcc = self.s * cc + self.w @ Ac + lm * mc
        if self.kernel == 'langmuir':
            flo, fhi = _f4_range(np.outer(p, blo), np.outer(p, bhi))
            g = (1.0 - self.rho)[:, None]; dAlo, dAhi = g * flo, g * fhi
            dmlo, dmhi = lm * blo / Hn, lm * bhi / Hn
            vlo, vhi = blo * r, bhi * r
            if self.kind == 'psi':
                dclo = (r / P) * (1.0 / (1.0 + vhi) - _ell(vlo)) + vlo / (1.0 + vlo)
                dchi = (r / P) * (1.0 / (1.0 + vlo) - _ell(vhi)) + vhi / (1.0 + vhi)
            else:
                f1, f2 = _f4_range(vlo, vhi); k = 1.0 - r / P
                dclo, dchi = (k * f1, k * f2) if k >= 0 else (k * f2, k * f1)
        else:
            act = (bc[None, :] > 1.0 / P) & (bc[None, :] < 1.0 / p[:, None])
            dAlo = np.where(act, np.outer(p, blo), 0.0); dAhi = np.where(act, np.outer(p, bhi), 0.0)
            on = bc > 1.0 / P
            dmlo, dmhi = np.where(on, lm * blo / Hn, 0.0), np.where(on, lm * bhi / Hn, 0.0)
            sat = bc * r > 1.0
            if self.kind == 'psi':
                g1, g2 = -np.log(np.maximum(blo * r, 1e-300)) / (blo * P), -np.log(np.maximum(bhi * r, 1e-300)) / (bhi * P)
                dclo = np.where(on, np.where(sat, 1.0, r * blo), np.where(sat, np.minimum(g1, g2), 0.0))
                dchi = np.where(on, np.where(sat, 1.0, r * bhi), np.where(sat, np.maximum(g1, g2), 0.0))
            else:
                dclo = np.where(sat, np.where(on, 0.0, -1.0 / (blo * P)), np.where(on, r * blo, 0.0))
                dchi = np.where(sat, np.where(on, 0.0, -1.0 / (bhi * P)), np.where(on, r * bhi, 0.0))
        if self.s < 0: dclo, dchi = -dchi, -dclo
        Dlo = dclo + self.wp @ dAlo + self.wn @ dAhi + dmlo
        Dhi = dchi + self.wp @ dAhi + self.wn @ dAlo + dmhi
        lbc = rcc - half * np.maximum(np.abs(Dlo), np.abs(Dhi)) - self._slack(np.abs(cc), self.m(bhi))
        return np.maximum(lbc, self.mono(blo, bhi)), rcc, bc

    def _far_lb(self, lo, hi, Bfar):
        return min(float(self.centered(lo, hi)[0].min()), self.tail(Bfar))

    def certify(self, thr, maxdepth=60, maxact=20000):
        """Return (v, cand, dlm).  After raising the mass multiplier by dlm >= 0 (any multiplier vector >= 0 is dual
        admissible; the dual objective drops by dlm), rc(b) >= -v for EVERY b in [0, inf).  cand lists affinities
        inside the column range with rc < -thr (new columns).  The region beyond the column range [B0, inf) is never
        given columns: it is certified by interval bounds on twelve further decades, an analytic tail bound, and,
        if needed, the smallest increase of the mass multiplier that makes that region dual feasible."""
        r, p = self.r, self.p
        kinks = np.array(list(1.0 / p) + [1.0 / r, math.e / r]) if self.kernel == 'hinge' else np.array([])
        B0 = float(self.S.bf[-1])
        near = np.unique(np.concatenate([self.S.bf[1:], kinks[kinks < B0]]))
        far = np.unique(np.concatenate([B0 * np.logspace(0, 12, 97), kinks[kinks > B0], [10.0 * max(1.0 / r, 1.0 / p[0], B0)]]))
        flo, fhi, Bfar = far[:-1], far[1:], float(far[-1])
        self.dlm = 0.0
        if self._far_lb(flo, fhi, Bfar) < -0.5 * thr:
            lm0 = self.lm; d = 1e-300
            while d < 1e300:                                   # exponential search, then bisection in log space
                self.lm = lm0 + d
                if self._far_lb(flo, fhi, Bfar) >= -0.5 * thr: break
                d *= 1e3
            a, b_ = d / 1e3, d
            for _ in range(40):
                mid = math.sqrt(a * b_); self.lm = lm0 + mid
                if self._far_lb(flo, fhi, Bfar) >= -0.5 * thr: b_ = mid
                else: a = mid
            self.dlm = b_; self.lm = lm0 + b_
        worst = min(self._far_lb(flo, fhi, Bfar), float(self.mono(np.array([0.0]), near[:1])[0]))
        lo, hi = near[:-1], near[1:]
        cb, cv = [], []
        for depth in range(maxdepth + 1):
            LB, RC, bc = self.centered(lo, hi)
            bad = RC < -thr                                    # true violations: candidates for new columns
            if bad.any(): cb.append(bc[bad]); cv.append(RC[bad])
            act = (LB < -thr) & ~bad                           # not yet certified, no violation found: bisect
            if (~act).any(): worst = min(worst, float(LB[~act].min()))
            if not act.any(): break
            if depth == maxdepth or int(act.sum()) > maxact:
                worst = min(worst, float(LB[act].min())); break
            mid = bc[act]
            lo, hi = np.concatenate([lo[act], mid]), np.concatenate([mid, hi[act]])
        self.depth = depth
        cand = np.array([])
        if cb:
            cb = np.concatenate(cb); cv = np.concatenate(cv); order = np.argsort(cv); pick = []
            for j in order:
                if all(abs(math.log(cb[j] / q)) > 0.02 for q in pick): pick.append(float(cb[j]))
                if len(pick) >= 12: break
            cand = np.array(pick)
        return max(0.0, -worst), cand, self.dlm


class StSource:
    def __init__(self, p, q, eps, H=None, kappa=None, per_decade=8, e=None, tol=1e-6, fine=16, kernel='langmuir'):
        self.kernel = kernel
        self.p = np.asarray(p, float); self.q = np.asarray(q, float)
        self.n = len(self.p); self.P = self.p[-1]; self.qm = self.q.max()
        self.eps = eps; self.tol = tol
        self.e = np.full(self.n, eps * self.qm) if e is None else np.asarray(e, float)
        self.Kmin = None; self.nlp = 0; self.ncol_added = 0; self.max_violation = 0.0
        self.max_depth = 0; self.n_tail_repair = 0; self.max_tail_repair = 0.0; self.n_unconverged = 0; self.n_fallback = 0
        if H is None:
            self.Kmin = min_henry(self.p, self.q, eps, kernel=kernel)
            if not math.isfinite(self.Kmin):
                self.feasible = False; self.H = float('nan'); return
            H = kappa * self.Kmin
        self.H = float(H)
        s_ref = max(self.q[0] / self.p[0], 1e-3 * self.qm / self.P)
        kap = max(self.H / s_ref, 1.0)
        self.b = _grid(self.p[0], self.P, kap, per_decade)
        self.bf = _grid(self.p[0], self.P, kap * 100, per_decade * fine, lo_dec=7.0)   # check grid
        if kernel == 'hinge':
            self.b = np.unique(np.concatenate([self.b, 1.0 / self.p])); self.bf = np.unique(np.concatenate([self.bf, 1.0 / self.p]))
        self.Hn = self.H / self.qm
        self.b_ub = np.concatenate([(self.q + self.e) / self.qm, -(self.q - self.e) / self.qm, [1.0]])
        self.Mu = float((self.q[-1] + self.e[-1]) / self.qm)      # bound on sum_j u_j from the last data row
        self.A_ub, self.dcol = self._cols(self.b)
        self.Af, self.dchk = self._cols(self.bf)
        res = linprog(np.zeros(len(self.b)), A_ub=self.A_ub, b_ub=self.b_ub, bounds=(0, None), method='highs')
        self.feasible = (res.status == 0)

    def _cols(self, b):
        """Columns in the variable u_j = contribution of site j to q(P_max)/q_max (well scaled for every Henry cap).
        Data rows: phi_b(p_i)/phi_b(P_max) in (0, 1]; mass row: 1/(Hn phi_b(P_max))."""
        ph = _phi(self.p, b, self.kernel)
        pn = ph[-1, :]
        A = ph / pn[None, :]
        return np.vstack([A, -A, (1.0 / (self.Hn * pn))[None, :]]), pn

    # ----- LP with column generation and a dual certificate that holds on the whole affinity half-line.
    # For multipliers lam >= 0 of the 2n+1 rows, weak duality gives, for EVERY feasible measure (not only grid ones),
    #   objective >= D + int rc(b) du(b) >= D - v * Mu,   D = -lam.b_ub,   v >= -inf_{b>=0} rc(b),   int du <= Mu,
    # where rc(b) = s*g_r(b)/phi_b(P) + sum_i (lam+_i - lam-_i) phi_b(p_i)/phi_b(P) + lam_m/(Hn phi_b(P)).
    # _certify() bounds inf_b rc(b) rigorously by interval arithmetic with adaptive bisection (plus both tails).
    maxit = 30          # column-generation rounds (1 = no column generation)
    repair = True       # subtract the dual repair term (False = raw grid LP value, an inner approximation)
    LPOPT = dict(primal_feasibility_tolerance=1e-9, dual_feasibility_tolerance=1e-9)
    certify = True      # True: continuum certificate (branch and bound); False: check-grid certificate only (legacy)

    def _obj(self, kind, r, s, b):
        """Objective per unit Henry mass on affinities b: s*psi_r(b) (kind 'psi') or s*phi_b(r) (kind 'phi')."""
        return s * (_kern(r, b, self.kernel) if kind == 'psi' else _qk(r, b, self.kernel))

    def _min(self, kind, r, s, maxit=None):
        maxit = self.maxit if maxit is None else maxit
        if not self.certify:
            cf = self._obj(kind, r, s, self.bf) / self.dchk
        res = None
        for it in range(maxit):
            c = self._obj(kind, r, s, self.b) / self.dcol
            # objective scale: the solver works with c/sig so that its absolute tolerances act as relative ones;
            # multipliers are mapped back, and the certificate below is evaluated in the original units.
            sig = float(abs(c[0]) if s < 0 else np.abs(c).min()); sig = sig if sig > 0 else 1.0
            res = linprog(c / sig, A_ub=self.A_ub, b_ub=self.b_ub, bounds=(0, None), method='highs', options=self.LPOPT)
            self.nlp += 1
            if res.status != 0:
                # numerical trouble on badly scaled records: fall back to looser settings.  This cannot hurt soundness,
                # because the certificate below is valid for ANY nonnegative multipliers, optimal or not.
                cs = 1.0 / np.maximum(1.0, np.abs(self.A_ub).max(axis=0))     # column equilibration (row duals are unchanged by it)
                for sg, sc, kw in ((sig, None, dict(method='highs')), (1.0, None, dict(method='highs')), (sig, cs, dict(method='highs')), (1.0, cs, dict(method='highs')),
                                   (1.0, None, dict(method='highs-ds')), (sig, cs, dict(method='highs-ds')), (1.0, None, dict(method='highs-ipm')), (sig, cs, dict(method='highs-ipm'))):
                    if sc is None: res = linprog(c / sg, A_ub=self.A_ub, b_ub=self.b_ub, bounds=(0, None), **kw)
                    else: res = linprog(c * sc / sg, A_ub=self.A_ub * sc[None, :], b_ub=self.b_ub, bounds=(0, None), **kw)
                    if res.status == 0: sig = sg; self.n_fallback += 1; break
                if res.status != 0:
                    raise RuntimeError('LP status %d' % res.status)
            lam = sig * np.maximum(-res.ineqlin.marginals, 0.0)   # multipliers of all 2n+1 rows, clipped at 0
            dual = float(-lam @ self.b_ub)                  # dual objective (not the primal value)
            ref = max(abs(dual), 1e-12)
            if self.certify:
                R = _RC(self, kind, r, s, lam)
                v, cand, dlm = R.certify(0.5 * self.tol * ref / self.Mu)
                dual -= dlm                                 # mass multiplier raised by dlm (right-hand side 1)
                self.max_depth = max(self.max_depth, R.depth); self.n_tail_repair += int(dlm > 0); self.max_tail_repair = max(self.max_tail_repair, dlm / ref)
            else:
                red = cf + self.Af.T @ lam                  # reduced cost on the check grid only
                v = float(max(0.0, -red.min()))
                neg = red < -0.1 * self.tol * ref / self.Mu
                loc = neg & (red <= np.roll(red, 1)) & (red <= np.roll(red, -1))
                idx = np.where(loc)[0]
                if len(idx) == 0: idx = np.array([int(np.argmin(red))])
                cand = self.bf[idx[np.argsort(red[idx])][:12]]
            gap = v * self.Mu                               # every feasible u has int du <= Mu
            if gap + (dlm if self.certify else 0.0) <= self.tol * ref or it == maxit - 1:
                break
            newb = np.setdiff1d(cand, self.b)
            if len(newb) == 0:
                break
            self.b = np.concatenate([self.b, newb])
            cols, pn = self._cols(newb)
            self.A_ub = np.hstack([self.A_ub, cols]); self.dcol = np.concatenate([self.dcol, pn])
            self.ncol_added += len(newb)
        tot = gap + (dlm if self.certify else 0.0)          # total certificate margin subtracted from the dual value
        self.max_violation = max(self.max_violation, tot / ref)
        self.last_violation = tot / ref
        self.n_unconverged += int(tot > self.tol * ref)
        if not self.repair:
            return float(res.fun) * sig, res
        return dual - gap, res

    RCAP = 1e4        # LPs are solved for r <= RCAP * Pmax; beyond that a sound closed-form continuation is used

    def _U_lp(self, r):
        val, _ = self._min('psi', r, -1.0)
        return -val * self.qm

    def _L_lp(self, r):
        val, _ = self._min('psi', r, 1.0)
        return val * self.qm

    def U(self, r):
        rc = self.RCAP * self.P
        if r <= rc:
            return self._U_lp(r)
        if not hasattr(self, '_Ucap'): self._Ucap = (self._U_lp(rc), self.q_hi(rc))
        return self._Ucap[0] + self._Ucap[1] * (r - rc) / rc      # q(t)/t is nonincreasing: q(t)/t <= q_hi(rc)/rc

    def L(self, r):
        rc = self.RCAP * self.P
        if r <= rc:
            return self._L_lp(r)
        if not hasattr(self, '_Lcap'): self._Lcap = (self._L_lp(rc), self.q_lo(rc))
        return self._Lcap[0] + self._Lcap[1] * math.log(r / rc)   # q nondecreasing

    def q_lo(self, r):
        val, _ = self._min('phi', r, 1.0)
        return val * self.qm

    def q_hi(self, r):
        val, _ = self._min('phi', r, -1.0)
        return -val * self.qm

    def extremal_sites(self, r, which='U', thresh=1e-7):
        """Sites (affinity, capacity, Henry mass) of an extremal measure for the envelope at r."""
        sgn = -1.0 if which == 'U' else 1.0
        _, res = self._min('psi', r, sgn)
        u = res.x; idx = np.where(u > thresh * max(u.sum(), 1e-300))[0]
        out = []
        for j in idx:
            hm = float(u[j] * self.qm / self.dcol[j])
            out.append((float(self.b[j]), hm / self.b[j] if self.b[j] > 0 else float('inf'), hm))
        return sorted(out)
