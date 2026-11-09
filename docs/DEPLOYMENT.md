# GitHub Pages deployment

The page is a standalone static website in `docs/`. It needs no npm install,
build service, API key or database.

For the repository owner, the first publication requires one repository setting:

1. Open **Settings → Pages**.
2. Under **Build and deployment → Source**, select **GitHub Actions**.
3. In **Actions → Project page**, run the workflow on `main` if the initial
   push preceded enablement.

Subsequent changes to `docs/` deploy automatically. The workflow publishes only
`docs/`, not the source package or local data. The project page is published at:

https://oliverdou776.github.io/What-Do-Tabulated-Isotherms-Certify-/

To preview locally:

```bash
python -m http.server 8000 --directory docs
```

Then open http://localhost:8000. The page uses local MP4 files, poster images
and optional WebVTT captions. The explorer under `docs/explorer/` performs table lookups in the browser; video playback performs no computation.
Video regeneration is documented in `media/README.md`.
