# Workbench recordings

These MP4s show actual interactions with the [public isotherm certifier](https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/explorer/). Each clip has an animated GIF, a poster, visible captions, click-target cursor overlays and an English WebVTT track. Screen captures are held and edited for readability; playback duration is not computation time.

| Recording | User workflow | Practical value and source |
|---|---|---|
| [Screen materials](screen-materials.mp4) | Ask whether CO₂/N₂ selectivity reaches 50 at 1 bar and 15% CO₂, run the query over 90 real pairs, filter the undetermined ones, inspect one material, loosen the Henry cap and export the CSV. | A shortlist that holds for every isotherm consistent with the tables, and the materials a best-fit screen would have passed without support. The browser performs lookups in envelope tables compiled by the Python package. |
| [Trust an extrapolation](trust-extrapolation.mp4) | Pick a real record whose upper 90% is hidden, compare six fits with the certified band, reveal the hidden measurements, and repeat on a record where the band misses some points. | See how far a single isotherm supports conclusions beyond its last point. Bands and fits were computed from the visible points only; records are a random sample of the paper's test split. |
| [Plan a measurement](plan-measurement.mp4) | Start from two zeolite 13X records cut at 1 bar and an undetermined target, read the lookahead advice, add the real CO₂ points and then the real N₂ points above 1 bar. | Learn which measurement decides a query before making it. The added points are the real ones from the source record; intervals after each extension were compiled by the Python package. |

Certified answers are relative to IAST, the isotherm class, the tolerance and the Henry cap. They are not predictions for real mixtures. `../explorer/data/reference.json` holds values from the Python engine that `tests/test_explorer.mjs` compares with the browser engine.

## Opening cartoon

`hero.mp4` / `hero.gif` is the sixteen-second looping summary at the top of the project page and the README: flue gas meets a zeolite that keeps CO₂ and lets N₂ pass, two isotherm tables are measured up to 1 bar, fitted models disagree, IsoCert certifies an interval and names the measurement that settles the question. The drawing is schematic and was authored as a CSS-animated SVG by Gemini 4 Argon from a written brief (`scripts/hero_story.html`); every number in it is computed from the zeolite 13X records shown in the planning view. Render it again with `python scripts/render_hero.py` (headless Chromium, frame-accurate seeking, FFmpeg).

## Repeat the workflow

From the repository root, serve the page with `python -m http.server 8000 --directory docs` and open `http://localhost:8000/explorer/`. Operate the controls listed above.

To rebuild the explorer data from ISODB records (about three CPU-hours; downloaded files are cached in the given directory):

```sh
python -m pip install -e '.[demo]'
python scripts/build_explorer_data.py --raw isodb_cache --only screen,extrapolate,measure,reference
python scripts/build_figures.py
```

To record and encode the clips again (FFmpeg with libx264 required):

```sh
python -m playwright install chromium
python scripts/record_demo.py --out captures
python scripts/render_demos.py --captures captures
```

`record_demo.py` drives the real page in a headless browser and stores, for each clip, the screenshots and a `take.json` with timestamps, captions and the positions of the operated elements. `render_demos.py` adds the caption band and the pointer overlay and writes H.264 MP4, an 840-pixel GIF, a PNG poster, WebVTT captions and a SHA-256 manifest.
