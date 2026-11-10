# Curated results

`published_results.csv` holds the manuscript's three result tables in long form: Table II (held-out extrapolation),
Table III (known truth) and Table IV (real pairs). `published_results.json` maps numbers quoted in the manuscript to
short keys (for example `HoCovAll` = 97.0%, the held-out coverage of the Langmuir-site class).

These are publication-level summaries that match the paper in `docs/assets/paper.pdf`. They are not raw run logs,
and they are separate from the explorer data in `docs/explorer/data/`, which are recomputed by
`scripts/build_explorer_data.py` for a sample of records.
