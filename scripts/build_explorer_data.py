"""Build the data behind the interactive explorer (docs/explorer/data/) from NIST ISODB records.

Every number shown by the explorer is produced here by the `isocert` package: compiled envelope tables, certified
loading bands, least-squares fits of six families, and lookahead advice.  The browser only performs table lookups and
IAST root finding on these tables (the paper's online stage).  The record selection is frozen in demo_sources.json.

    python scripts/build_explorer_data.py --raw DIR      # DIR caches ISODB JSON files; missing ones are downloaded
"""
import argparse
import json
import math
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from isocert import Record, interval, next_measurement
from isocert.baselines import fit_all
from isocert.compiled import _down, _up
from isocert.isodb import fetch_isotherm, load_isotherm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs' / 'explorer' / 'data'
CAPS = (10.0, 100.0, 1000.0)
GAS = {'Carbon Dioxide': 'CO2', 'Nitrogen': 'N2', 'Methane': 'CH4', 'Ethene': 'C2H4', 'Ethane': 'C2H6', 'Hydrogen': 'H2',
       'Oxygen': 'O2', 'Acetylene': 'C2H2', 'Propene': 'C3H6', 'N-propane': 'C3H8', 'Argon': 'Ar', 'Krypton': 'Kr',
       'Xenon': 'Xe', 'Carbon monoxide': 'CO'}
RAW = None


def sig(v, n=5):
    return float('%.*g' % (n, v))


def load(name):
    path = RAW / (name.replace('/', '_') + '.json')
    if path.exists():
        return load_isotherm(path)
    rec = fetch_isotherm(name)                          # cache the upstream file next to the others
    return rec


def fits_of(p, q, eps):
    """Six least-squares fits: parameters, full-scale RMSE, max residual, and whether the fit is within tolerance."""
    out = {}
    for name, (m, rmse, sup) in fit_all(p, q).items():
        out[name] = dict(theta=[float('%.9g' % t) for t in m.theta], rmse=sig(rmse, 4), sup=sig(sup, 4), ok=bool(sup <= eps))
    return out


def band_of(R, lo, hi, n=36):
    grid = np.geomspace(lo, hi, n)
    ql, qh = R.band(grid)
    return dict(p=[sig(x, 6) for x in grid], lo=[_down(float(x)) for x in ql], hi=[_up(float(x)) for x in qh])


def table_of(R):
    d = R.compile().to_dict()
    return dict(U=d['U'], L=d['L'], q=d['q_lo'], H=d['H'], qt=d['q_hi_top']), d['r']


def meta_of(rec):
    return dict(file=rec['name'], gas=GAS.get(rec['gas'], rec['gas']), adsorbent=rec['adsorbent'], T=rec['T'], doi=rec['doi'])


# ------------------------------------------------------------------ view 1: screen a set of materials
def screen_record(name):
    try:
        rec = load(name); R = Record(rec['p'], rec['q'])
        out = dict(meta_of(rec), p=[sig(x, 6) for x in R.p], q=[sig(x, 6) for x in R.q], eps=sig(R.eps, 4), floor=sig(R.eps_star, 4),
                   fits=fits_of(R.p, R.q, R.eps), band=band_of(R, R.p[0] / 3, R.p[-1] * 30), tables={})
        for cap in CAPS:
            t, r = table_of(R.with_cap(cap)); out['tables']['%g' % cap] = t; out['r'] = r
        return name, out
    except Exception as ex:
        return name, dict(error=repr(ex)[:200])


# ------------------------------------------------------------------ view 2: trust an extrapolation
def extrapolate_record(name):
    try:
        rec = load(name); p, q = rec['p'], rec['q']; vis = p <= p[-1] / 10.0
        R = Record(p[vis], q[vis]); H = Record(p[vis], q[vis], kernel='hinge', eps=R.eps)
        e = R.eps * float(R.q.max()); ph, qh = p[~vis], q[~vis]
        lo, hi = R.band(ph); hlo, hhi = H.band(ph)
        return name, dict(meta_of(rec), p=[sig(x, 6) for x in R.p], q=[sig(x, 6) for x in R.q], eps=sig(R.eps, 4), e=sig(e, 5),
                          hidden=dict(p=[sig(x, 6) for x in ph], q=[sig(x, 6) for x in qh], lo=[_down(float(x)) for x in lo], hi=[_up(float(x)) for x in hi],
                                      hlo=[_down(float(x)) for x in hlo], hhi=[_up(float(x)) for x in hhi]),
                          band=band_of(R, R.p[0], p[-1] * 1.05, 44), hinge=band_of(H, R.p[0], p[-1] * 1.05, 44), fits=fits_of(R.p, R.q, R.eps))
    except Exception as ex:
        return name, dict(error=repr(ex)[:200])


# ------------------------------------------------------------------ view 3: plan a measurement
def levels(p, n=6):
    """Indices of the last point of each extension level: measured up to 1 bar, then real points above it."""
    base = int(np.sum(p <= 1.0)); extra = len(p) - base
    steps = sorted(set(int(round(x)) for x in np.linspace(0, extra, min(n, extra + 1))))
    return [base + s for s in steps]


def measure_level(args):
    name, n_pts, y, P = args
    try:
        rec = load(name); R = Record(rec['p'][:n_pts], rec['q'][:n_pts])
        t, r = table_of(R)
        return (name, n_pts), dict(n=n_pts, p_max=sig(float(rec['p'][n_pts - 1]), 5), eps=sig(R.eps, 4), floor=sig(R.eps_star, 4), admitted=R.admitted,
                                   table=t, r=r, fits=fits_of(R.p, R.q, R.eps), band=band_of(R, rec['p'][0] / 3, rec['p'][-1] * 3))
    except Exception as ex:
        return (name, n_pts), dict(n=n_pts, error=repr(ex)[:200])


def measure_advice(args):
    a, b, y, P = args
    try:
        ra, rb = load(a), load(b)
        A = Record(ra['p'][ra['p'] <= 1.0], ra['q'][ra['p'] <= 1.0]); B = Record(rb['p'][rb['p'] <= 1.0], rb['q'][rb['p'] <= 1.0])
        base = interval(A, B, P, y)
        adv = next_measurement(A, B, P, y)
        return (a, b), dict(base=[sig(base.lnS_lo, 6), sig(base.lnS_hi, 6)],
                            advice=[dict(gas=d['gas'], p=sig(d['pressure'], 4), lo=_down(math.log(d['S_lo'])) if d['S_lo'] > 0 else None,
                                         hi=_up(math.log(d['S_hi'])) if math.isfinite(d['S_hi']) else None, w=sig(d['width'], 5)) for d in adv])
    except Exception as ex:
        return (a, b), dict(error=repr(ex)[:200])


def write_reference():
    """Values computed by the Python engine from the stored tables and fits; tests/test_explorer.mjs compares the browser engine with them."""
    from isocert import CompiledRecord, decide
    from isocert.baselines import FAMILIES
    from isocert.iast import point_selectivity
    fam = {c.name: c for c in FAMILIES}
    tab = lambda r, t: CompiledRecord(r, t['U'], t['L'], t['q'], t['H'], t['qt'])
    best = lambda fits: (lambda k: fam[k](fits[k]['theta']))(min(fits, key=lambda k: fits[k]['rmse']))
    cases = []
    for slug in ('co2-n2', 'co2-ch4', 'c2h4-c2h6'):
        d = json.loads((OUT / ('screen-%s.json' % slug)).read_text())
        for pr in d['pairs'][:6]:
            for P, y, cap in ((d['P'], d['y'], '10'), (0.2, 0.5, '100'), (8.0, 0.05, '1000'), (1e-5, 0.3, '10'), (2e4, 0.5, '10')):
                A, B = tab(pr['A']['r'], pr['A']['tables'][cap]), tab(pr['B']['r'], pr['B']['tables'][cap])
                iv = interval(A, B, P, y); tau = d['tau']
                fin = lambda v: v if math.isfinite(v) else ('inf' if v > 0 else '-inf')
                cases.append(dict(view='screen', file=slug, id=pr['id'], P=P, y=y, cap=cap, tau=tau, lo=fin(iv.lnS_lo), hi=fin(iv.lnS_hi), answer=decide(A, B, P, y, tau),
                                  best=math.log(point_selectivity(best(pr['A']['fits']), best(pr['B']['fits']), P, y))))
    d = json.loads((OUT / 'measure.json').read_text())
    for pr in d['pairs'][:5]:
        for i, j in ((0, 0), (len(pr['A']['levels']) - 1, len(pr['B']['levels']) - 1)):
            LA, LB = pr['A']['levels'][i], pr['B']['levels'][j]
            iv = interval(tab(LA['r'], LA['table']), tab(LB['r'], LB['table']), pr['P'], pr['y'])
            cases.append(dict(view='measure', id=pr['id'], i=i, j=j, lo=iv.lnS_lo, hi=iv.lnS_hi))
    (OUT / 'reference.json').write_text(json.dumps(dict(source='Python engine (isocert) evaluated on the stored tables and fits', cases=cases), indent=0))
    print('reference cases', len(cases))


def main():
    global RAW
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--raw', type=Path, required=True, help='directory with ISODB isotherm JSON files (name with "/" replaced by "_")')
    ap.add_argument('--procs', type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument('--only', default='screen,extrapolate,measure')
    args = ap.parse_args(); RAW = args.raw
    src = json.loads((ROOT / 'scripts' / 'demo_sources.json').read_text()); OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    defaults = {'Carbon Dioxide/Nitrogen': ('co2-n2', 0.15, 50.0), 'Carbon Dioxide/Methane': ('co2-ch4', 0.5, 5.0), 'Ethene/Ethane': ('c2h4-c2h6', 0.5, 2.0)}
    with Pool(args.procs) as pool:
        if 'screen' in args.only:
            names = sorted(set(f for v in src['screen'].values() for pr in v for f in pr))
            done = dict(pool.imap_unordered(screen_record, names, chunksize=1)); print('screen records', len(done), round(time.time() - t0), flush=True)
            for pair, lst in src['screen'].items():
                slug, y, tau = defaults[pair]
                pairs = [dict(id=i + 1, A=done[a], B=done[b]) for i, (a, b) in enumerate(lst) if 'error' not in done[a] and 'error' not in done[b]]
                gA, gB = [GAS[g] for g in pair.split('/')]
                (OUT / ('screen-%s.json' % slug)).write_text(json.dumps(dict(gasA=gA, gasB=gB, y=y, P=1.0, tau=tau, caps=[int(c) for c in CAPS], pairs=pairs), separators=(',', ':')))
                print(pair, len(pairs), 'pairs; errors', [done[f]['error'] for pr in lst for f in pr if 'error' in done[f]][:2], flush=True)
        if 'extrapolate' in args.only:
            done = dict(pool.imap_unordered(extrapolate_record, src['extrapolate'], chunksize=1))
            recs = [dict(id=i + 1, **done[n]) for i, n in enumerate(src['extrapolate']) if 'error' not in done[n]]
            (OUT / 'extrapolate.json').write_text(json.dumps(dict(records=recs), separators=(',', ':'))); print('extrapolate', len(recs), round(time.time() - t0), flush=True)
        if 'measure' in args.only:
            tasks, advs, plan = [], [], []
            for a, b in src['measure']:
                ra, rb = load(a), load(b); y = 0.15 if rb['gas'] == 'Nitrogen' else 0.5
                la, lb = levels(ra['p']), levels(rb['p']); plan.append((a, b, y, la, lb, ra, rb))
                tasks += [(a, n, y, 1.0) for n in la] + [(b, n, y, 1.0) for n in lb]; advs.append((a, b, y, 1.0))
            lev = dict(pool.imap_unordered(measure_level, tasks, chunksize=1)); print('measure levels', len(lev), round(time.time() - t0), flush=True)
            adv = dict(pool.imap_unordered(measure_advice, advs, chunksize=1))
            out = []
            for i, (a, b, y, la, lb, ra, rb) in enumerate(plan):
                if any('error' in lev[(a, n)] for n in la) or any('error' in lev[(b, n)] for n in lb) or 'error' in adv[(a, b)]:
                    print('skip', a, b, [lev[(a, n)].get('error') for n in la], [lev[(b, n)].get('error') for n in lb], adv[(a, b)].get('error')); continue
                def side(rec, ls, name):
                    return dict(meta_of(rec), p=[sig(x, 6) for x in rec['p']], q=[sig(x, 6) for x in rec['q']], levels=[lev[(name, n)] for n in ls])
                out.append(dict(id=i + 1, y=y, P=1.0, tau=50.0 if y == 0.15 else 5.0, A=side(ra, la, a), B=side(rb, lb, b), **adv[(a, b)]))
            (OUT / 'measure.json').write_text(json.dumps(dict(pairs=out), separators=(',', ':'))); print('measure', len(out), round(time.time() - t0), flush=True)
    if 'reference' in args.only:
        write_reference()
    for f in sorted(OUT.glob('*.json')): print(f.name, f.stat().st_size // 1024, 'KB')


if __name__ == '__main__':
    main()
