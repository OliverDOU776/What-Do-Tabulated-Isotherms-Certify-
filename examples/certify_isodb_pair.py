"""Certify a selectivity target for two NIST ISODB records (same adsorbent and temperature).

    python examples/certify_isodb_pair.py 10.1021Ie2014982.Isotherm13 10.1021Ie2014982.Isotherm14 --P 1 --y 0.15 --target 50

Arguments are ISODB filenames (downloaded from the public API) or paths to local ISODB JSON files.
"""
import argparse
import os

from isocert import Record, decide, interval
from isocert.baselines import fit_all
from isocert.iast import point_selectivity
from isocert.isodb import fetch_isotherm, load_isotherm


def get(arg):
    rec = load_isotherm(arg) if os.path.exists(arg) else fetch_isotherm(arg)
    if rec is None: raise SystemExit(f'{arg}: not a pure-gas record in bar with a supported loading unit')
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('record_a'); ap.add_argument('record_b')
    ap.add_argument('--P', type=float, default=1.0, help='total pressure (bar)')
    ap.add_argument('--y', type=float, default=0.15, help='gas-phase mole fraction of gas A')
    ap.add_argument('--target', type=float, default=None, help='selectivity threshold')
    ap.add_argument('--kappa', type=float, default=10.0, help='Henry-cap multiplier')
    a = ap.parse_args()
    ra, rb = get(a.record_a), get(a.record_b)
    A = Record(ra['p'], ra['q'], kappa=a.kappa, name=ra['gas']); B = Record(rb['p'], rb['q'], kappa=a.kappa, name=rb['gas'])
    print(f"{ra['adsorbent']}  ({ra['T']:g} K, doi:{ra['doi']})")
    for R in (A, B):
        print(f'  {R.name}: {len(R.p)} points up to {R.p[-1]:g} bar, misfit floor {100 * R.eps_star:.2f}%, admitted {R.admitted}')
    if not (A.admitted and B.admitted): print('  warning: a record lies outside the Langmuir-site class at the 3% gate')
    fa, fb = fit_all(A.p, A.q), fit_all(B.p, B.q)
    ba, bb = min(fa, key=lambda k: fa[k][1]), min(fb, key=lambda k: fb[k][1])
    iv = interval(A, B, a.P, a.y)
    print(f'Best fits ({ba}, {bb}): S = {point_selectivity(fa[ba][0], fb[bb][0], a.P, a.y):.2f}')
    print(f'Certified interval at P = {a.P:g} bar, y = {a.y:g}: S in [{iv.S_lo:.2f}, {iv.S_hi:.2f}]')
    print(f'IAST reads {B.name} at up to {max(iv.r_B):.3g} bar (last measured point {B.p[-1]:g} bar)')
    if a.target: print(f'S >= {a.target:g}?  {decide(A, B, a.P, a.y, a.target)}')


if __name__ == '__main__':
    main()
