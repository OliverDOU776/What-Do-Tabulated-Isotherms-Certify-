<div align="center">

# What Do Tabulated Isotherms Certify?<br>Sharp Certain Answers to Mixture-Adsorption Screening Queries

<p><a href="mailto:zijia.wang@said.oxford.edu">Zijia Wang</a><sup>∗</sup>, <a href="mailto:haoliang.hong@eng.ox.ac.uk">Haoliang Hong</a><sup>∗</sup>, <a href="mailto:lixuhui25@mails.ucas.ac.cn">Xuhui Li</a><sup>†</sup>, <a href="mailto:j.bai22@imperial.ac.uk">Jingyang Bai</a><sup>‡</sup>, <a href="mailto:binjian.nie@eng.ox.ac.uk">Binjian Nie</a><sup>∗</sup></p>

<sup>∗</sup>University of Oxford, United Kingdom · <sup>†</sup>University of Chinese Academy of Sciences, China · <sup>‡</sup>Imperial College London, United Kingdom

<details><summary>Author contact</summary>

Zijia Wang: [zijia.wang@said.oxford.edu](mailto:zijia.wang@said.oxford.edu)

Haoliang Hong: [haoliang.hong@eng.ox.ac.uk](mailto:haoliang.hong@eng.ox.ac.uk)

Xuhui Li: [lixuhui25@mails.ucas.ac.cn](mailto:lixuhui25@mails.ucas.ac.cn)

Jingyang Bai: [j.bai22@imperial.ac.uk](mailto:j.bai22@imperial.ac.uk)

Binjian Nie: [binjian.nie@eng.ox.ac.uk](mailto:binjian.nie@eng.ox.ac.uk)

</details>

[Paper](docs/assets/paper.pdf) · [Project page](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/) · [Interactive demo](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/explorer/) · [Demo video](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/media/screen-materials.mp4) · [BibTeX](#citation)

**Find out what two isotherm tables certify about a gas mixture, and which measurement would settle the rest.**

</div>

![Sixteen-second cartoon: flue gas meets a zeolite that keeps CO2 and lets N2 pass; two isotherm tables measured up to 1 bar; models that fit both give selectivities of 217 to 2,018; IsoCert certifies 178 to 15,497 (undetermined for a target of 300) and advises measuring N2 above 1 bar; with those points the interval is 421 to 2,810, a certain yes](docs/media/hero.gif?v=20261109-2)

<sub>Schematic drawing; every number is computed from real records.</sub>

Take CO₂ and N₂ on zeolite 13X at 298 K, both measured up to 1 bar ([source record](https://doi.org/10.1016/j.micromeso.2015.05.044)). Is the flue-gas selectivity at least 300? **Twenty-four pairs of models fit both tables within tolerance and answer 217 to 2,018.** Mixture theory reads the N₂ isotherm far beyond its last measured point, where those fits disagree. IsoCert computes what holds for *every* isotherm consistent with the tables: 178 to 15,497, so the answer is **UNDETERMINED**. It advises measuring N₂ above 1 bar. With the N₂ points that the source measured up to 10 bar, the interval becomes 421 to 2,810 and the answer is a certain **YES**. Adding the CO₂ points up to 10 bar instead changes nothing.

![The same example in data: N2 isotherm on zeolite 13X measured to 1 bar with six fits and the certified band; fitted models give selectivities of 217 to 2,018, the tables certify 178 to 15,497, and after adding N2 points to 10 bar the interval is 421 to 2,810](docs/assets/opening-example.png)

On 607 NIST ISODB records measured up to 10 bar, certified bands computed from a low-pressure prefix of each record cover **97.0%** of 12,402 held-out measurements and 96.5% at 5–10× extrapolation, where the best-fit model falls within ±3% of only 6.9%. Across 1,930 real gas pairs, at targets a factor 2 below or above the best-fit selectivity, the tables certify 66% of pairs as YES and 50% as NO, and lookahead advice settles 68% of undetermined pairs in a median of 2 added measurements. The guarantee is relative to IAST, an isotherm class, a tolerance and a Henry cap, and the paper reports where it stops holding.

## See how to use it

[![Screen recording of running a selectivity query over real records, filtering the undetermined materials, inspecting one and exporting the certain answers](docs/media/screen-materials.gif?v=20261109-1)](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/media/screen-materials.mp4)

**Ask a screening question, read the certain answers, then settle what is open.** The recording operates the same [public workbench](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/explorer/) that you can open yourself. It runs on real records.

1. Choose the gas pair, total pressure, gas composition and the selectivity you need.
2. Run the query. Every material is YES, NO or UNDETERMINED for all isotherms consistent with its two tables.
3. Select a material to see its tables, the certified bands, the fitted models and where IAST reads each record.
4. Export the answers as CSV, or take an undetermined material to the measurement planner.

In the recorded query (CO₂/N₂, 1 bar, 15% CO₂, *S* ≥ 50), a best-fit screen shortlists 33 of 90 real pairs; **14 of those are not certain**, and 19 pairs certainly meet the target.

| View | What you do | What you learn |
|---|---|---|
| [Screen materials](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/explorer/#screen) · [Video](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/media/screen-materials.mp4) | Set mixture, pressure, composition, target and Henry cap; inspect and export | Which adsorbents certainly pass, which certainly fail, and which a fit would pass without support |
| [Trust an extrapolation](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/explorer/#extrapolate) · [Video](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/media/trust-extrapolation.mp4) | Compare fits with the certified band on a record whose upper 90% is hidden, then reveal it | How far one isotherm supports conclusions beyond its last point, including a record where the band misses |
| [Plan a measurement](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/explorer/#measure) · [Video](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/media/plan-measurement.mp4) | Read the advice for an undetermined query; add the real points above 1 bar for either gas | Which gas and pressure range decides the query, before measuring |

The workbench evaluates the paper's online stage in the browser: lookups in envelope tables that the Python package compiled from 411 ISODB records, and binary IAST on those tables. To certify your own tables, follow [Installation](#installation) and [Usage](#usage).

## Method

A fit replaces a table by one curve. IsoCert keeps every curve the table allows.

A record with points (pᵢ, qᵢ), tolerance eᵢ = ε·maxⱼ qⱼ and Henry slope cap H = κ·K_min defines a cone of site mixtures q(p) = c₀p + ∫ φ_b(p) dμ(b), with c₀ ≥ 0, a nonnegative site measure μ, data boxes qᵢ − eᵢ ≤ q(pᵢ) ≤ qᵢ + eᵢ and q′(0) ≤ H. Langmuir sites φ_b(p) = p/(1 + bp) give the Langmuir-site Stieltjes class C_L, which contains single- and multi-site Langmuir, Tóth and UNILAN isotherms with finite Henry constants; hinge sites min(p, 1/b) give the shape-only class C_h of all nonnegative, nondecreasing, concave isotherms. The dual-site Langmuir set C₂ ⊆ C_L ⊂ C_h brackets them from inside. The spreading pressure is linear in the site measure, so its upper and lower envelopes U(r), L(r) over a class are linear semi-infinite programs. The **Envelope Reduction Theorem** uses the monotonicity of the IAST root in the two spreading pressures: the exact selectivity interval over both classes is

**S⁺ = IAST(U_A, L_B),  S⁻ = IAST(L_A, U_B),**

both ends are attained, and each envelope value is attained by an isotherm with at most n + 2 sites. A threshold query S ≥ τ is YES if S⁻ ≥ τ, NO if S⁺ < τ, and UNDETERMINED otherwise; it is decided from four envelope values.

Two laws bound what a table can certify. **Extrapolation-Horizon Law:** relative uncertainty stays of order ε inside the measured window but becomes of order one at the extrapolation ratio ρ = r/P_max ∼ (1 + 1/c)/ε = Θ(1/ε), where c is the saturation degree of a Langmuir isotherm at P_max. **Log-Henry Law:** below the lowest measured pressure p₁, the spreading-pressure width grows as Θ(q₁⁻ ln H), so an explicit Henry cap is necessary; each decision reports its critical Henry cap κ⋆, the largest multiplier under which it still holds.

Offline, IsoCert admits a record when its misfit floor ε⋆, the optimum of a linear program, is at most 3%, solves each envelope by column generation with a dual certificate that is valid for every site affinity, and compiles the record once into monotone outer-bound tables. Online, it answers a threshold query from four table lookups and an interval query by a scalar root solve on the tables. For an undetermined query, a one-step lookahead ranks candidate measurements by the interval width they are predicted to leave.

## Experimental results

**Admission.** The misfit-floor LP admits 68.5% of the 13,183 ISODB isotherms (9,033 records, 14 gases) into the Langmuir-site class at ε⋆ ≤ 3%, with a median misfit floor of 1.54%.

**Held-out extrapolation on real records.** Methods see only a low-pressure prefix of each record, up to a cut-off pressure P_cut. Coverage (%) of the hidden points and median band ratio q_hi/q_lo, by extrapolation ratio ρ = p/P_cut.

| Method | ρ ∈ (1, 2] | ρ ∈ (2, 5] | ρ ∈ (5, 10] | all ρ ∈ (1, 10] | band ratio (all) |
|---|---:|---:|---:|---:|---:|
| Best-fit model (point) | 62.3 | 19.5 | 6.9 | 19.6 | 1.00 |
| Six-family ensemble | 80.1 | 61.6 | 58.8 | 63.0 | 1.24 |
| Dual-site parameter set | 97.4 | 95.6 | 94.5 | 95.3 | 2.43 |
| Hinge class (shape only) | 99.0 | 98.6 | 99.1 | 98.9 | 4.29 |
| **Langmuir-site class (ours)** | **98.3** | **97.3** | **96.5** | **97.0** | **2.51** |

607 Type-I records with P_max ≤ 10 bar, 12,402 hidden points. A point counts as covered inside the band widened by the record's tolerance, ε = 3% of the maximum loading; across 412 groups of replicate measurements (1,184 curves), the median inter-replicate discrepancy is 1.42% of the maximum loading and this tolerance covers 97.1% of repeat curves. At 5–10× extrapolation the Langmuir-site band is 1.9× tighter than the shape-only band (median ratio 3.68 against 6.85).

**Real pairs.** Across 1,930 pairs and 12 gas separations, the certified interval of ln S has median width 1.25 (a factor of 3.5). A six-family ensemble spans 26% of that width, and its spread correlates with fit error at 0.21 (Spearman); residual-bootstrap intervals lie strictly inside the certified interval on 51% of pairs. The set of all dual-site Langmuir parameters within tolerance spans a median 96%: the ambiguity is already present in ordinary two-site models, and point fitting hides it. At thresholds a factor 2 below and above the best-fit selectivity, IsoCert certifies 66% of pairs as YES and 50% as NO.

**Known truth.** On 990 synthetic instances, the interval contains the true selectivity in all 670 that respect the Henry cap, with 0 false decisions among 1,592 certified threshold answers.

**Measurement advice.** On 100 undetermined pairs, one-step lookahead reaches the target precision for 68% in a median of 2 added measurements and cuts the median width by 94%; adding points inside the measured window reaches it for 7%.

**Cost.** Compiling a pair of records takes 9.9 s and pays off after 8 queries. Afterwards an interval query takes 0.15 ms (8118× faster than direct LP evaluation at 1235 ms) and a threshold query 25 µs (2437× faster).

**Where the guarantee ends.** Answers are relative to IAST, the class, the tolerance and the cap. On all 272 binary mixture series in ISODB (2,512 measured points), the certified interval contains the measured selectivity for 60% at κ = 10 (75% at κ = 1000; a fit ensemble: 24%). For records measured above 10 bar, held-out coverage drops to 84.2%. Over all 990 synthetic instances, including those whose Henry constant violates the cap, coverage is 95.1% and 1.8% of calls at a 1.25× margin (1.3% at 2×) are false. The dual-site parameter set reaches a median 96% of the certified width; what the certified class adds is the outer guarantee. The measurement advice is a myopic one-step heuristic, not an optimal experimental design.

[Result tables II–IV](evidence/published_results.csv) · [Numbers quoted in the paper](evidence/published_results.json) · [Full paper](docs/assets/paper.pdf)

## Installation

The implementation requires **Python 3.10+, NumPy and SciPy** (HiGHS linear programming). No GPU is required.

```bash
git clone https://github.com/OliverDOU776/What-Do-Tabulated-Isotherms-Certify-.git
cd What-Do-Tabulated-Isotherms-Certify-
python -m venv .venv
source .venv/bin/activate
python -m pip install .
python examples/quickstart.py
```

## Usage

```python
from isocert import Record, interval, decide, critical_cap, next_measurement

co2 = Record(p_co2, q_co2)                 # pressures (bar) and loadings (mmol/g) of one table
n2 = Record(p_n2, q_n2)

iv = interval(co2, n2, P=1.0, y=0.15)      # certified [S-, S+] at 1 bar, 15% CO2
print(iv.S_lo, iv.S_hi)
print(decide(co2, n2, P=1.0, y=0.15, tau=50))            # 'YES' | 'NO' | 'UNDETERMINED'
print(critical_cap(co2, n2, P=1.0, y=0.15, tau=50))      # largest Henry-cap multiplier keeping the answer
print(next_measurement(co2, n2, P=1.0, y=0.15, top=3))   # which gas, which pressure
```

`Record` computes the misfit floor of the table and admits it if the floor is at most 3% of full scale (`record.admitted`). By default it sets the tolerance to max(3 × misfit floor, 1%) of full scale and the Henry cap to κ = 10 times the smallest Henry constant consistent with the table; `eps=0.03` gives the fixed 3% tolerance stated in the paper, and `eps`, `kappa`, `H` and `kernel='hinge'` override the defaults. Any consistent pressure and loading units work. `record.band(pressures)` returns the certified loading band.

Compile once, then query by lookup; a compiled table is a small JSON file that needs no solver:

```python
table_a, table_b = co2.compile(), n2.compile()
table_a.save('co2.table.json')
print(interval(table_a, table_b, P=0.5, y=0.1))
```

```bash
python examples/certify_isodb_pair.py RECORD_A RECORD_B --P 1 --y 0.15 --target 50   # two NIST ISODB filenames
python examples/compile_and_screen.py RECORDS_DIR --gas-a "Carbon Dioxide" --gas-b Nitrogen --target 50
```

`isocert.baselines` provides what current practice reports for comparison: least-squares fits of six families and the dual-site parameter set.

## Reproducibility

```bash
python -m pip install '.[test]'
python -m pytest -q
node --test tests/test_explorer.mjs
```

The Python tests check, on isotherms whose site affinities lie on no solver grid, that envelopes contain the generating isotherm, that intervals contain the true selectivity and decisions never contradict it, that the extrapolation-horizon witness bounds the lower envelope, that classes are nested, and that compiled tables are outer bounds that survive a JSON round trip. The Node test checks that the browser engine reproduces the Python engine on the stored tables. The numerical kernel is the one used for the paper's experiments; it reproduces the paper's stored intervals exactly on sampled pairs.

The explorer data, the opening figure and the recordings are regenerated by the scripts described in [docs/media/README.md](docs/media/README.md). The repository provides functional code and compact derived results, without internal run logs or research administration records. The explorer covers a sample of records; it does not reproduce the complete paper campaign.

| Directory | Contents |
|---|---|
| `src/isocert/` | Envelope solver with continuum certificate, compiled tables, IAST on envelopes, record interface, baselines, ISODB reader |
| `examples/` | Quickstart, certifying two ISODB records, compiling and screening a folder |
| `scripts/` | Explorer data build, opening cartoon and figure, UI recording and encoding |
| `tests/` | Functional checks of the Python engine and of the browser engine |
| `docs/` | GitHub Pages source: project page, interactive workbench, paper, figures and videos |
| `evidence/` | Curated publication-level result tables |

The tests have finite coverage. Certificates hold up to double-precision rounding; verified interval arithmetic is not used. Certified answers do not establish experimental IAST accuracy or physical validity for real mixtures. See the paper for the full model and evaluation scope.

## Citation

```bibtex
@misc{wang2026isotherms,
  title  = {What Do Tabulated Isotherms Certify? Sharp Certain Answers to Mixture-Adsorption Screening Queries},
  author = {Wang, Zijia and Hong, Haoliang and Li, Xuhui and Bai, Jingyang and Nie, Binjian},
  year   = {2026},
  url    = {https://github.com/OliverDOU776/What-Do-Tabulated-Isotherms-Certify-}
}
```

## Attribution

Software and manuscript: Zijia Wang, Haoliang Hong, Xuhui Li, Jingyang Bai, and Binjian Nie (contacts above). Isotherm records shown in the explorer come from the NIST/ARPA-E adsorption database and retain their upstream terms; see [NOTICE.md](NOTICE.md). The project page is an original implementation informed by the organization of [VGGT](https://vgg-t.github.io/), [ReAct](https://react-lm.github.io/), [WebArena](https://webarena.dev/og/) and [Self-Refine](https://selfrefine.info/).
