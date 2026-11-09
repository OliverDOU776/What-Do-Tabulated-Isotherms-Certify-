"""Render the opening cartoon (scripts/hero_story.html, CSS-animated SVG) to docs/media/hero.{mp4,gif} and a poster.

The page is opened in headless Chromium; every CSS animation is paused and moved to the same time stamp before each
screenshot, so the export is deterministic.

    python scripts/render_hero.py            # needs playwright (chromium) and FFmpeg with libx264
"""
import argparse
import subprocess
import tempfile
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
MEDIA = ROOT / 'docs' / 'media'
SEEK = '(t) => { for (const a of document.getAnimations()) { a.pause(); a.currentTime = t; } }'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--source', type=Path, default=ROOT / 'scripts' / 'hero_story.html')
    ap.add_argument('--duration', type=float, default=16.0); ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--poster', type=float, default=15.2, help='time (s) of the poster frame')
    ap.add_argument('--stills', type=Path, help='also save one still per second into this directory (for checking)')
    a = ap.parse_args()
    n = round(a.duration * a.fps)
    with tempfile.TemporaryDirectory() as tmp, sync_playwright() as p:
        browser = p.chromium.launch(); page = browser.new_page(viewport={'width': 1280, 'height': 640}, device_scale_factor=1)
        page.goto(a.source.resolve().as_uri()); page.wait_for_timeout(300)
        for i in range(n):
            page.evaluate(SEEK, 1000.0 * i / a.fps); page.screenshot(path=f'{tmp}/{i:04d}.png')
        page.evaluate(SEEK, 1000.0 * a.poster); page.screenshot(path=str(MEDIA / 'hero-poster.png'))
        if a.stills:
            a.stills.mkdir(parents=True, exist_ok=True)
            for s in range(int(a.duration)):
                page.evaluate(SEEK, 1000.0 * (s + 0.9)); page.screenshot(path=str(a.stills / f'{s:02d}.png'))
        browser.close()
        ff = ['ffmpeg', '-y', '-v', 'error', '-framerate', str(a.fps), '-i', f'{tmp}/%04d.png']
        subprocess.run(ff + ['-an', '-c:v', 'libx264', '-preset', 'medium', '-crf', '18', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', str(MEDIA / 'hero.mp4')], check=True)
        flt = 'fps=10,scale=840:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=96[p];[b][p]paletteuse=dither=bayer:bayer_scale=4'
        subprocess.run(ff + ['-filter_complex', flt, '-loop', '0', str(MEDIA / 'hero.gif')], check=True)
    for f in ('hero.mp4', 'hero.gif', 'hero-poster.png'): print(f, (MEDIA / f).stat().st_size // 1024, 'KB')


if __name__ == '__main__':
    main()
