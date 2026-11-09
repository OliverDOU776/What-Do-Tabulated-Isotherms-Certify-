"""Certify a threshold query on two tabulated isotherms (synthetic numbers; replace them with your own tables).

    python examples/quickstart.py
"""
import numpy as np
from isocert import Record, CompiledRecord, interval, decide, critical_cap, next_measurement
from isocert.baselines import fit_all
from isocert.iast import point_selectivity

# Two pure-gas records on one adsorbent, measured up to 1 bar (pressure in bar, loading in mmol/g).
p = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0])
q_co2 = np.array([0.36, 0.63, 1.17, 1.67, 2.20, 2.61, 2.86, 3.07, 3.19, 3.28])
q_n2 = np.array([0.004, 0.008, 0.020, 0.039, 0.077, 0.131, 0.183, 0.249, 0.297, 0.343])

co2, n2 = Record(p, q_co2, name='CO2'), Record(p, q_n2, name='N2')
for r in (co2, n2):
    print(f'{r.name}: misfit floor {100 * r.eps_star:.2f}%  tolerance {100 * r.eps:.1f}%  admitted {r.admitted}')

P, y = 1.0, 0.15                                    # flue gas: 1 bar, 15% CO2
fits_a, fits_b = fit_all(co2.p, co2.q), fit_all(n2.p, n2.q)
print('\nWhat fitted models report (all within tolerance on both tables):')
for name in fits_a:
    if name in fits_b and fits_a[name][2] <= co2.eps and fits_b[name][2] <= n2.eps:
        print(f'  {name:<10s} S = {point_selectivity(fits_a[name][0], fits_b[name][0], P, y):7.1f}')

iv = interval(co2, n2, P, y)
print(f'\nCertified interval: S in [{iv.S_lo:.1f}, {iv.S_hi:.1f}]')
print(f'IAST reads the N2 record at {min(iv.r_B):.1f}-{max(iv.r_B):.1f} bar; its last measured point is {n2.p[-1]:g} bar.')
for tau in (50.0, 80.0):                            # every fit above says "yes" to both targets
    print(f'  S >= {tau:g}?  {decide(co2, n2, P, y, tau)}')
print('Critical Henry-cap multiplier of the answer to S >= 50:', critical_cap(co2, n2, P, y, 50.0, kappa_max=1e3, steps=4))

print('\nWhich single measurement would narrow the interval most?')
for a in next_measurement(co2, n2, P, y, candidates=[0.001, 0.003, 0.3, 3.0, 10.0, 30.0], top=3):
    gas = 'CO2' if a['gas'] == 'A' else 'N2'
    print(f"  {gas} at {a['pressure']:g} bar  ->  predicted [{a['S_lo']:.1f}, {a['S_hi']:.1f}]")

# Compile once, then answer any later query by table lookup (no solver needed).
table_a, table_b = co2.compile(), n2.compile()
table_a.save('co2_table.json'); table_a = CompiledRecord.load('co2_table.json')
for P2, y2 in ((1.0, 0.15), (0.5, 0.15), (2.0, 0.05)):
    c = interval(table_a, table_b, P2, y2)
    print(f'lookup  P = {P2:g} bar, y = {y2:g}:  S in [{c.S_lo:.1f}, {c.S_hi:.1f}]')
