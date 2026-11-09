"""Read isotherm files of the NIST/ARPA-E adsorption database (ISODB, https://adsorption.nist.gov/isodb/)."""
import json
import math
import urllib.parse
import urllib.request
import numpy as np

API = 'https://adsorption.nist.gov/isodb/api/isotherm/'
MOLAR = {'Carbon Dioxide': 44.01, 'Nitrogen': 28.014, 'Methane': 16.043, 'Ethane': 30.07, 'Ethene': 28.054,
         'N-propane': 44.097, 'Propene': 42.081, 'Oxygen': 31.998, 'Argon': 39.948, 'Krypton': 83.798,
         'Xenon': 131.293, 'Carbon monoxide': 28.01, 'Acetylene': 26.038, 'Hydrogen': 2.016}


def to_mmol_g(unit, gas):
    """Factor converting a loading unit to mmol/g, or None when the unit is not supported."""
    u = unit.strip().lower().replace(' ', '')
    if u in ('mmol/g', 'mol/kg'): return 1.0
    if u in ('mmol/kg', 'umol/g', 'µmol/g'): return 1e-3
    if u in ('cm3(stp)/g', 'ml(stp)/g', 'cc(stp)/g', 'ml/g', 'cm3/g'): return 1.0 / 22.414
    if u == 'mg/g' and gas in MOLAR: return 1.0 / MOLAR[gas]
    if u == 'g/g' and gas in MOLAR: return 1000.0 / MOLAR[gas]
    return None


def parse_isotherm(d):
    """Pure-gas ISODB JSON (dict) -> dict(name, gas, adsorbent, doi, T, p [bar], q [mmol/g]); None if unsupported."""
    if len(d.get('adsorbates', [])) != 1 or str(d.get('pressureUnits', '')).lower() != 'bar': return None
    gas = d['adsorbates'][0].get('name', '')
    conv = to_mmol_g(d.get('adsorptionUnits', ''), gas)
    if conv is None: return None
    agg = {}
    for x in d.get('isotherm_data', []):
        try: p = float(x['pressure']); q = float(x['species_data'][0]['adsorption'])
        except (KeyError, IndexError, TypeError, ValueError): continue
        if p > 0 and math.isfinite(q): agg.setdefault(p, []).append(q)
    if len(agg) < 3: return None
    p = np.array(sorted(agg)); q = np.array([np.mean(agg[x]) for x in p]) * conv
    return dict(name=d.get('filename', ''), gas=gas, adsorbent=d.get('adsorbent', {}).get('name', ''), doi=d.get('DOI', ''),
                T=float(d.get('temperature', 'nan')), p=p, q=q)


def load_isotherm(path):
    with open(path) as f: return parse_isotherm(json.load(f))


def fetch_isotherm(filename, timeout=40):
    """Download one isotherm by its ISODB filename.  Data remain subject to the database's and the source article's terms."""
    req = urllib.request.Request(API + urllib.parse.quote(filename) + '.json', headers={'User-Agent': 'isocert/0.1'})
    with urllib.request.urlopen(req, timeout=timeout) as f: return parse_isotherm(json.loads(f.read()))
