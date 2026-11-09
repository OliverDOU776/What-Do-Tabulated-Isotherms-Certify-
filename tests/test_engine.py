"""Functional checks of the certified envelopes, intervals, compiled tables and the high-level interface."""
import json
import math

import numpy as np
import pytest

from isocert import CompiledRecord, Record, critical_cap, decide, interval, next_measurement
from isocert.baselines import DualSiteSet, fit_all
from isocert.iast import point_selectivity, solve_z


class DualSite:
    """Ground-truth dual-site Langmuir isotherm with affinities that lie on no solver grid."""
    def __init__(self, m, b): self.m, self.b = np.asarray(m, float), np.asarray(b, float)
    def q(self, p): return float((self.m * p / (1 + self.b * p)).sum())
    def Pi(self, r): return float((self.m * np.log1p(self.b * r) / self.b).sum())
    @property
    def henry(self): return float(self.m.sum())


P_A = np.geomspace(0.01, 1.0, 12)
P_B = np.geomspace(0.02, 1.0, 10)
TRUE_A = DualSite([37.3, 2.9], [11.7, 0.313])
TRUE_B = DualSite([0.83, 0.21], [0.71, 0.0457])


@pytest.fixture(scope='module')
def pair():
    A = Record(P_A, [TRUE_A.q(p) for p in P_A], eps=0.01, H=1.5 * TRUE_A.henry, name='A')
    B = Record(P_B, [TRUE_B.q(p) for p in P_B], eps=0.01, H=1.5 * TRUE_B.henry, name='B')
    return A, B


def test_envelopes_contain_the_generating_isotherm(pair):
    for R, truth in zip(pair, (TRUE_A, TRUE_B)):
        for r in (1e-3, 0.03, 0.5, 2.0, 40.0, 3000.0):
            assert R.L(r) <= truth.Pi(r) <= R.U(r)
            lo, hi = R.band([r])
            assert lo[0] <= truth.q(r) <= hi[0]
        assert R.envelope.max_violation < 1e-5          # certificate margin subtracted from every bound


def test_interval_contains_true_selectivity_and_decisions_agree(pair):
    A, B = pair
    for P, y in ((1.0, 0.15), (0.2, 0.5), (5.0, 0.3)):
        s_true = point_selectivity(TRUE_A, TRUE_B, P, y)
        iv = interval(A, B, P, y)
        assert iv.S_lo <= s_true <= iv.S_hi
        assert decide(A, B, P, y, iv.S_lo * 0.98) == 'YES'
        assert decide(A, B, P, y, iv.S_hi * 1.02) == 'NO'
        assert decide(A, B, P, y, s_true) == 'UNDETERMINED'
        assert iv.answer(iv.S_lo * 0.98) == 'YES' and iv.answer(s_true) == 'UNDETERMINED'


def test_extrapolation_horizon_witness():
    """Theorem 3: for a Henry-like table the one-site witness bounds L from above, for any number of points."""
    K, Pmax, eps, a = 2.0, 1.0, 0.02, 0.0
    b_a = (a + eps) / ((1 - eps) * Pmax)
    for n in (6, 24):
        p = np.linspace(Pmax / n, Pmax, n)
        R = Record(p, K * p, eps=eps, H=K)
        for rho in (3.0, 30.0, 300.0):
            r = rho * Pmax
            witness = (1 + a) * K * math.log1p(b_a * r) / b_a
            assert R.L(r) <= witness * (1 + 1e-9)
            assert R.L(r) >= 0.95 * witness              # and the LP value is close to it
            assert R.U(r) >= K * r * (1 - 1e-9)


def test_hinge_class_contains_langmuir_class(pair):
    A, B = pair
    hA = Record(A.p, A.q, kernel='hinge', eps=A.eps, H=A.H); hB = Record(B.p, B.q, kernel='hinge', eps=B.eps, H=B.H)
    lang, hinge = interval(A, B, 1.0, 0.15), interval(hA, hB, 1.0, 0.15)
    assert hinge.lnS_lo <= lang.lnS_lo + 2e-3 and hinge.lnS_hi >= lang.lnS_hi - 2e-3


def test_compiled_tables_are_outer_bounds_and_survive_json(pair, tmp_path):
    A, B = pair
    CA, CB = A.compile(), B.compile()
    CA.save(tmp_path / 'a.json'); CA2 = CompiledRecord.load(tmp_path / 'a.json')
    CB2 = CompiledRecord.from_dict(json.loads(json.dumps(CB.to_dict())))
    for P, y in ((1.0, 0.15), (0.05, 0.5), (20.0, 0.8), (1e-4, 0.5), (500.0, 0.2)):
        d, c = interval(A, B, P, y), interval(CA2, CB2, P, y)
        assert c.lnS_lo <= d.lnS_lo + 1e-9 and c.lnS_hi >= d.lnS_hi - 1e-9
    d, c = interval(A, B, 1.0, 0.15), interval(CA2, CB2, 1.0, 0.15)
    assert c.width - d.width < 0.05                       # tables add little width inside the tabulated range
    assert decide(CA2, CB2, 1.0, 0.15, d.S_lo * 0.9) == 'YES'


def test_looser_cap_widens_and_critical_cap_is_bracketed(pair):
    A = Record(P_A, [TRUE_A.q(p) for p in P_A], eps=0.01); B = Record(P_B, [TRUE_B.q(p) for p in P_B], eps=0.01)
    i10, i100 = interval(A, B, 1.0, 0.15), interval(A.with_cap(100), B.with_cap(100), 1.0, 0.15)
    assert i100.lnS_lo <= i10.lnS_lo + 2e-3 and i100.lnS_hi >= i10.lnS_hi - 2e-3
    tau = i10.S_lo * 0.7
    k = critical_cap(A, B, 1.0, 0.15, tau, kappa_max=1e3, steps=4)
    assert k >= 10.0 and decide(A.with_cap(k), B.with_cap(k), 1.0, 0.15, tau) == 'YES'
    assert critical_cap(A, B, 1.0, 0.15, math.sqrt(i10.S_lo * i10.S_hi)) is None


def test_new_measurement_cannot_widen_and_advice_is_sorted(pair):
    A, B = pair
    base = interval(A, B, 1.0, 0.15)
    B2 = B.with_point(8.0, TRUE_B.q(8.0))
    after = interval(A, B2, 1.0, 0.15)
    assert after.lnS_lo >= base.lnS_lo - 2e-3 and after.lnS_hi <= base.lnS_hi + 2e-3 and after.width < base.width
    adv = next_measurement(A, B, 1.0, 0.15, candidates=[0.003, 0.3, 3.0, 30.0])
    assert [a['width'] for a in adv] == sorted(a['width'] for a in adv)
    assert adv[0]['width'] <= base.width + 2e-3


def test_baselines_fit_and_dual_site_set_is_inside_the_certified_class(pair):
    A, B = pair
    fits = fit_all(A.p, A.q)
    assert fits['DSL'][2] < 1e-3 and set(fits) >= {'Langmuir', 'DSL', 'Sips', 'Toth', 'UNILAN', 'Freundlich'}
    DA, DB = DualSiteSet(A.p, A.q, A.eps, A.H), DualSiteSet(B.p, B.q, B.eps, B.H)
    assert not DA.empty and not DB.empty
    cert = interval(A, B, 1.0, 0.15)
    z_lo = solve_z(DA.L, DB.U, 1.0, 0.15, side='lo'); z_hi = solve_z(DA.U, DB.L, 1.0, 0.15, side='up'); c = math.log(0.15 / 0.85)
    assert cert.lnS_lo - 2e-3 <= z_lo - c and z_hi - c <= cert.lnS_hi + 2e-3


def test_admission_and_invalid_input():
    p = np.geomspace(0.01, 1, 10)
    s_shaped = (3 * p) ** 3 / (1 + (3 * p) ** 3)               # not a site mixture: convex at low pressure
    assert not Record(p, s_shaped).admitted
    assert Record(p, 2 * p / (1 + 3 * p)).admitted
    with pytest.raises(ValueError): Record([0.1, 0.1, 0.2], [1, 2, 3])
    with pytest.raises(ValueError): Record([0.1, 0.2, 0.3], [0, 0, 0])
    with pytest.raises(ValueError): Record([0.1, 0.2], [1, 2])
