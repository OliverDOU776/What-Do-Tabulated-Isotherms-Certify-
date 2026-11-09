"""Binary IAST on envelopes: certified selectivity interval and the threshold predicate.

IAST equates the reduced spreading pressures of the two gases,
    Pi_A(P y / x) = Pi_B(P (1 - y) / (1 - x)),      S = (x / (1 - x)) / (y / (1 - y)),
where x is the adsorbed mole fraction of gas A.  Because the root is monotone in the two spreading pressures, the
extreme selectivities over two site-mixture classes are obtained from their envelopes:
    S+ = IAST(U_A, L_B),    S- = IAST(L_A, U_B).
Any objects exposing U(r) and L(r) can be passed (a Record, a CompiledRecord, or a parametric model via `as_envelope`).
"""
import math
from typing import NamedTuple


def solve_z(Pi1, Pi2, Ptot, y1, tol=1e-3, zlo=-30.0, zhi=30.0, side=None):
    """Root of Pi1(P y1 (1+e^-z)) - Pi2(P y2 (1+e^z)), which is decreasing in z = ln(x1/x2).
    side='up' returns an upper bound of the root (right end of the final bracket, +inf if the root lies beyond the
    search range), side='lo' a lower bound (left end, -inf beyond the range); side=None returns the midpoint."""
    y2 = 1.0 - y1
    f = lambda z: Pi1(Ptot * y1 * (1.0 + math.exp(-z))) - Pi2(Ptot * y2 * (1.0 + math.exp(z)))
    if f(zlo) <= 0:
        return -math.inf if side == 'lo' else zlo
    if f(zhi) >= 0:
        return math.inf if side == 'up' else zhi
    while zhi - zlo > tol:
        mid = 0.5 * (zlo + zhi)
        if f(mid) > 0: zlo = mid
        else: zhi = mid
    if side == 'up': return zhi
    if side == 'lo': return zlo
    return 0.5 * (zlo + zhi)


def lnS_bounds(A, B, Ptot, y1, tol=1e-3):
    """(ln S-, ln S+, z-, z+) from envelope objects A, B; z = ln(x/(1-x)) at the two interval ends."""
    c = math.log(y1 / (1.0 - y1))
    zu = solve_z(A.U, B.L, Ptot, y1, tol, side='up')
    zl = solve_z(A.L, B.U, Ptot, y1, tol, side='lo')
    return zl - c, zu - c, zl, zu


class Interval(NamedTuple):
    """Certified interval of the IAST selectivity S_AB at (P, y)."""
    S_lo: float
    S_hi: float
    lnS_lo: float
    lnS_hi: float
    r_A: tuple      # IAST pressures at which gas A is evaluated at the two ends (low end, high end)
    r_B: tuple

    @property
    def width(self):
        return self.lnS_hi - self.lnS_lo

    def answer(self, tau):
        if self.S_lo >= tau: return 'YES'
        if self.S_hi < tau: return 'NO'
        return 'UNDETERMINED'


def interval(A, B, P, y, tol=1e-3):
    """Certified interval [S-, S+] of the selectivity of gas A over gas B at total pressure P and gas-phase mole
    fraction y of A.  Every pair of isotherms in the two classes has its IAST selectivity inside the interval."""
    lo, up, zl, zu = lnS_bounds(_env(A), _env(B), P, y, tol)
    ex = lambda z: math.exp(z) if abs(z) < 700 else (math.inf if z > 0 else 0.0)
    return Interval(ex(lo), ex(up), lo, up, (P * y * (1 + ex(-zl)), P * y * (1 + ex(-zu))),
                    (P * (1 - y) * (1 + ex(zl)), P * (1 - y) * (1 + ex(zu))))


def decide(A, B, P, y, tau):
    """Certain answer to the threshold query "S_AB >= tau?": 'YES', 'NO' or 'UNDETERMINED'.
    Threshold predicate: with x_tau = tau*y/(1-y+tau*y), r_A = P*y/x_tau and r_B = P*(1-y)/(1-x_tau),
    YES iff L_A(r_A) >= U_B(r_B) and NO iff U_A(r_A) < L_B(r_B).  Four envelope values, no root solve."""
    A, B = _env(A), _env(B)
    z = math.log(tau) + math.log(y / (1.0 - y))
    rA = P * y * (1.0 + math.exp(-z)); rB = P * (1.0 - y) * (1.0 + math.exp(z))
    if A.L(rA) >= B.U(rB): return 'YES'
    if A.U(rA) < B.L(rB): return 'NO'
    return 'UNDETERMINED'


def point_selectivity(model_A, model_B, P, y):
    """IAST selectivity of two single isotherms given by objects with a Pi(r) method (e.g. parametric fits)."""
    z = solve_z(model_A.Pi, model_B.Pi, P, y, 1e-6, -60.0, 60.0)
    return math.exp(z - math.log(y / (1.0 - y)))


def _env(obj):
    return obj.envelope if hasattr(obj, 'envelope') else obj
