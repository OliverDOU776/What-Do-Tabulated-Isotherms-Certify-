"""Compile a folder of ISODB records once, then screen every same-adsorbent pair of two gases from the tables.

    python examples/compile_and_screen.py RECORDS_DIR --gas-a "Carbon Dioxide" --gas-b Nitrogen --P 1 --y 0.15 --target 50

Compiled tables are cached next to the records (*.table.json); later runs and other conditions only do lookups.
"""
import argparse
import collections
from pathlib import Path

from isocert import CompiledRecord, Record, decide, interval
from isocert.isodb import load_isotherm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('directory', type=Path)
    ap.add_argument('--gas-a', default='Carbon Dioxide'); ap.add_argument('--gas-b', default='Nitrogen')
    ap.add_argument('--P', type=float, default=1.0); ap.add_argument('--y', type=float, default=0.15); ap.add_argument('--target', type=float, default=50.0)
    a = ap.parse_args()
    groups = collections.defaultdict(dict)
    for path in sorted(a.directory.glob('*.json')):
        if path.name.endswith('.table.json'): continue
        rec = load_isotherm(path)
        if rec is None or rec['gas'] not in (a.gas_a, a.gas_b) or len(rec['p']) < 6: continue
        cache = path.with_suffix('.table.json')
        if cache.exists():
            table = CompiledRecord.load(cache)
        else:
            R = Record(rec['p'], rec['q'], name=rec['name'])
            if not R.admitted: continue                      # outside the class at the 3% gate
            table = R.compile(); table.save(cache)
        key = (rec['adsorbent'], rec['doi'], round(rec['T']))
        groups[key].setdefault(rec['gas'], table)
    print(f'{"adsorbent":<32s} {"T/K":>5s}  {"S-":>9s} {"S+":>9s}  answer to S >= {a.target:g}')
    for (ads, doi, T), g in sorted(groups.items()):
        if a.gas_a in g and a.gas_b in g:
            iv = interval(g[a.gas_a], g[a.gas_b], a.P, a.y)
            print(f'{ads[:32]:<32s} {T:>5d}  {iv.S_lo:>9.2f} {iv.S_hi:>9.2f}  {decide(g[a.gas_a], g[a.gas_b], a.P, a.y, a.target)}')


if __name__ == '__main__':
    main()
