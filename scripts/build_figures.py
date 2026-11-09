"""Draw the opening example of the project page from the explorer data (real ISODB records, no hand-set numbers).

    python scripts/build_figures.py        # writes docs/assets/opening-example.svg and prints the numbers quoted on the page
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from isocert import CompiledRecord, Record, interval
from isocert.baselines import FAMILIES
from isocert.iast import point_selectivity

ROOT = Path(__file__).resolve().parents[1]
PAIR, TARGET = 3, 300.0          # Zeolite 13X, CO2/N2, 298 K (doi:10.1016/j.micromeso.2015.05.044)
C = dict(ink='#30363d', soft='#67727d', rule='#dfe4e8', band='#7fb0d8', line='#2973ac', fit='#c0564b', read='#e2a03f', grey='#9aa5ae', green='#187f79')
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 10.5, 'svg.fonttype': 'none',
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': C['soft'], 'axes.labelcolor': C['ink'], 'xtick.color': C['soft'], 'ytick.color': C['soft']})


def main():
    fam = {c.name: c for c in FAMILIES}
    pr = next(p for p in json.loads((ROOT / 'docs/explorer/data/measure.json').read_text())['pairs'] if p['id'] == PAIR)
    tab = lambda L: CompiledRecord(L['r'], L['table']['U'], L['table']['L'], L['table']['q'], L['table']['H'], L['table']['qt'])
    A, B = pr['A'], pr['B']; A0, B0, B1 = A['levels'][0], B['levels'][0], B['levels'][-1]
    base = interval(tab(A0), tab(B0), pr['P'], pr['y']); after = interval(tab(A0), tab(B1), pr['P'], pr['y'])
    fits = lambda L: {k: fam[k](v['theta']) for k, v in L['fits'].items() if v['ok']}
    fa, fb = fits(A0), fits(B0)
    ens = sorted(point_selectivity(a, b, pr['P'], pr['y']) for a in fa.values() for b in fb.values())

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.4, 3.55), gridspec_kw=dict(width_ratios=[1, 1.2], wspace=0.16, left=0.06, right=0.985, top=0.87, bottom=0.165))
    # (a) the weaker gas: measured to 1 bar, read by IAST far beyond.  The band is recomputed here on a longer pressure range.
    n0 = B0['n']; p, q = np.array(B['p']), np.array(B['q']); rb = sorted(base.r_B); xmax = 4000.0
    R = Record(p[:n0], q[:n0]); grid = np.geomspace(p[0] * 0.7, xmax, 70); lo, hi = R.band(grid)
    ax.fill_between(grid, lo, hi, color=C['band'], alpha=0.38, lw=0); ax.plot(grid, lo, color=C['band'], lw=1); ax.plot(grid, hi, color=C['band'], lw=1)
    for m in fb.values(): ax.plot(grid, m.q(grid), color=C['fit'], lw=1.2)
    ax.axvspan(rb[0], min(rb[1], xmax), color=C['read'], alpha=0.18, lw=0)
    ax.plot(p[:n0], q[:n0], 'o', ms=5, color=C['ink'], mec='white', mew=0.6, zorder=5)
    ax.axvline(p[n0 - 1], color=C['soft'], ls=(0, (2, 3)), lw=1)
    ymax = 9.0; ax.set_xscale('log'); ax.set_xlim(p[0] * 0.7, xmax); ax.set_ylim(0, ymax)
    ax.set_xlabel('pressure (bar)'); ax.set_ylabel('N$_2$ loading (mmol/g)')
    ax.set_title('N$_2$ on zeolite 13X, measured up to 1 bar', loc='left', fontsize=11.5, color=C['ink'], pad=8)
    ax.text(p[n0 - 1] * 0.82, ymax * 0.90, 'last measured\npoint', ha='right', va='top', fontsize=9, color=C['soft'])
    ax.text(math.sqrt(rb[0] * min(rb[1], xmax)), ymax * 0.955, 'IAST reads the\nrecord here', ha='center', va='top', fontsize=9, color='#8a5d12')
    ax.text(0.035, 0.52, '%d fits within\ntolerance' % len(fb), transform=ax.transAxes, fontsize=9, color=C['fit'])
    ax.text(0.035, 0.36, 'certified band', transform=ax.transAxes, fontsize=9, color=C['line'])
    # (b) what each approach says about the selectivity
    rows = [('%d pairs of fitted models, all within tolerance on both tables' % len(ens), None), ('Certified from the tables (both gases measured up to 1 bar)', base),
            ('Certified after adding the N$_2$ points measured up to 10 bar', after)]
    x0 = 95.0
    for i, (lab, iv) in enumerate(rows):
        y = 2 - i; bx.text(x0, y + 0.33, lab, fontsize=10, color=C['ink'], va='bottom', bbox=dict(fc='white', ec='none', pad=1.5), zorder=4)
        if iv is None:
            bx.vlines(ens, y - 0.17, y + 0.17, color=C['fit'], lw=1.6)
            bx.text(ens[-1] * 1.3, y, 'S = {:,.0f} to {:,.0f}'.format(ens[0], ens[-1]), va='center', fontsize=10, color=C['fit'])
        else:
            bx.plot([iv.S_lo, iv.S_hi], [y, y], color=C['grey'] if i == 1 else C['line'], lw=9, solid_capstyle='butt')
            ans = iv.answer(TARGET)
            bx.text(iv.S_hi * 1.3, y, 'S = {:,.0f} to {:,.0f}'.format(iv.S_lo, iv.S_hi), va='center', fontsize=10, color=C['ink'])
            bx.text(0.985, y, ans, transform=bx.get_yaxis_transform(), ha='right', va='center', fontsize=10, fontweight='bold', color=C['green'] if ans == 'YES' else '#7d5c10')
    bx.axvline(TARGET, color=C['ink'], ls=(0, (4, 3)), lw=1, ymax=0.9); bx.text(TARGET * 1.08, -0.42, 'target: S ≥ %g' % TARGET, ha='left', fontsize=9.5, color=C['ink'])
    bx.set_xscale('log'); bx.set_xlim(85, 1.5e7); bx.set_ylim(-0.6, 2.9); bx.set_yticks([]); bx.spines['left'].set_visible(False)
    bx.set_xlabel('CO$_2$/N$_2$ IAST selectivity at 1 bar, 15% CO$_2$'); bx.set_title('Is the selectivity at least 300?', loc='left', fontsize=11.5, color=C['ink'], pad=8)
    out = ROOT / 'docs/assets/opening-example.svg'; fig.savefig(out); fig.savefig(out.with_suffix('.png'), dpi=170)
    print('fits within tolerance: %d pairs, S = %.0f-%.0f' % (len(ens), ens[0], ens[-1]))
    print('certified, data to 1 bar: [%.0f, %.0f] -> %s' % (base.S_lo, base.S_hi, base.answer(TARGET)))
    print('certified, N2 to %.3g bar: [%.0f, %.0f] -> %s' % (B1['p_max'], after.S_lo, after.S_hi, after.answer(TARGET)))
    print('IAST reads N2 at %.1f-%.0f bar; last measured point %.2f bar' % (rb[0], rb[1], p[n0 - 1]))


if __name__ == '__main__':
    main()
