"""Drive the real explorer in a headless browser and save timestamped captures for scripts/render_demos.py.

    python -m pip install playwright && python -m playwright install chromium
    python scripts/record_demo.py --out captures      # then: python scripts/render_demos.py --captures captures

Every frame is a screenshot of docs/explorer/ after a real interaction; cursor positions are the centres of the
elements that were operated.  Playback time is editing time, not computation time.
"""
import argparse
import functools
import http.server
import json
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
VIEW = {'width': 1200, 'height': 1090}


class Take:
    def __init__(self, page, directory):
        self.page, self.dir, self.t, self.frames, self.cues = page, directory, 0.0, [], []
        directory.mkdir(parents=True, exist_ok=True)

    def _shot(self, cursor=None):
        name = '%03d.png' % len(self.frames)
        self.page.screenshot(path=str(self.dir / name))
        frame = {'t': round(self.t, 2), 'file': name}
        if cursor: frame['cursor'] = cursor
        self.frames.append(frame)

    def say(self, caption, hold=4.0):
        self.cues.append({'t': round(self.t, 2), 'caption': caption}); self._shot(self.frames[-1].get('cursor') if self.frames else None); self.t += hold

    def _centre(self, selector, fx=0.5):
        el = self.page.locator(selector).first; el.scroll_into_view_if_needed(); box = el.bounding_box()
        return {'x': round(box['x'] + box['width'] * fx), 'y': round(box['y'] + box['height'] / 2)}

    def act(self, selector, action, caption=None, hold=4.0, fx=0.5, settle=350):
        """Show the pointer on the element, perform the real interaction, then capture the result."""
        c = self._centre(selector, fx); self._shot(c); self.t += 0.8
        action(self.page.locator(selector).first); self.page.wait_for_timeout(settle)
        c = dict(self._centre(selector, fx), click=True)
        if caption: self.cues.append({'t': round(self.t, 2), 'caption': caption})
        self._shot(c); self.t += hold

    def mark(self):
        self.poster = round(self.t - 1.5, 2)       # poster frame: the state just captured

    def save(self, poster_time=None):
        poster_time = getattr(self, 'poster', poster_time)
        (self.dir / 'take.json').write_text(json.dumps({'duration': round(self.t, 2), 'poster_time': poster_time, 'frames': self.frames, 'cues': self.cues}, indent=1))


def slide(value):
    return lambda el: el.evaluate('(e, v) => { e.value = v; e.dispatchEvent(new Event("input", {bubbles: true})); }', str(value))


def screen_materials(page, out, base):
    page.goto(base + '/explorer/#screen'); page.wait_for_selector('#rows tr'); page.evaluate('window.scrollTo(0, 62)')
    t = Take(page, out / 'screen-materials')
    t.say('A screening question: which adsorbents reach a CO₂/N₂ selectivity of 50 at 1 bar, 15% CO₂?', 4.5)
    t.act('#query button[type=submit]', lambda e: e.click(), 'Each row is a real pair of records. 19 of 90 certainly meet the target; 38 cannot be decided.', 5)
    t.act('#filters button[data-f=UNDETERMINED]', lambda e: e.click(), 'A best-fit screen would shortlist 33 materials. 14 of them are not certain.', 4.5)
    t.act('#rows tr:has-text("TbL") button', lambda e: e.click(), 'Inspect one. Every fitted model says yes, yet isotherms consistent with the same tables fall below 50.', 5.5)
    t.mark()
    t.say('The reason: IAST reads the N₂ record far beyond its last measured point (dashed lines).', 4.5)
    t.act('#filters button[data-f=YES]', lambda e: e.click(), 'These 19 answers hold for every isotherm consistent with the data.', 4)
    t.act('select[name=cap]', lambda e: e.select_option('1000'), None, 1.2)
    t.act('#query button[type=submit]', lambda e: e.click(), 'Loosen the Henry cap 100-fold: the answers that remain do not depend on it.', 5)
    with page.expect_download():
        t.act('#export', lambda e: e.click(), 'Export the certain answers together with the query conditions.', 4)
    t.save(9.5)


def trust_extrapolation(page, out, base):
    page.goto(base + '/explorer/#extrapolate'); page.wait_for_selector('#ex-plot svg'); page.evaluate('window.scrollTo(0, 62)')
    t = Take(page, out / 'trust-extrapolation')
    t.act('#ex-record', lambda e: e.select_option(label='PF-600 ZnCl2 · CO₂ · 273 K'), 'A real CO₂ record. Only its low-pressure tenth is shown; six models fit it within tolerance.', 5)
    t.say('Beyond the window the fits fan out. The blue band is every isotherm the visible points allow.', 4.5)
    t.act('#ex-reveal', lambda e: e.click(), 'Reveal the measurements: all 31 lie in the certified band; the best fit holds for 1.', 5.5)
    t.mark()
    t.act('#ex-hinge', lambda e: e.check(), 'Without the site structure, a shape-only band is also safe but much wider.', 4.5)
    t.act('#ex-record', lambda e: e.select_option(label='UiO-66 · N₂ · 298 K'), 'Another record, hidden again.', 3)
    t.act('#ex-reveal', lambda e: e.click(), 'Here 14 of 136 points escape the band: the guarantee is relative to the class, and the page says so.', 5.5)
    t.save(10.5)


def plan_measurement(page, out, base):
    page.goto(base + '/explorer/#measure'); page.wait_for_selector('#me-strip svg'); page.evaluate('window.scrollTo(0, 62)')
    t = Take(page, out / 'plan-measurement')
    t.act('#me-pair', lambda e: e.select_option(label='Zeolite 13X · CO₂/N₂ · 298 K'), 'Zeolite 13X, both records measured up to 1 bar. Is the CO₂/N₂ selectivity at least 300?', 4.5)
    t.act('#me-tau', lambda e: e.fill('300'), 'Well-fitting models answer 217 to 2,018. The tables alone cannot decide: UNDETERMINED.', 5.5)
    t.say('The advice: measure N₂ above 1 bar. More CO₂ points are predicted to change nothing.', 5)
    page.evaluate('window.scrollTo(0, 62)')
    n_a = int(page.locator('#me-a').get_attribute('max')); n_b = int(page.locator('#me-b').get_attribute('max'))
    t.act('#me-a', slide(n_a), 'Add the real CO₂ measurements up to 10 bar: the interval does not move.', 4.5, fx=0.97)
    t.act('#me-a', slide(0), None, 1.0, fx=0.03)
    t.act('#me-b', slide(max(1, n_b // 2)), 'Now add the real N₂ measurements, a few points at a time.', 3.5, fx=0.5)
    t.act('#me-b', slide(n_b), 'With N₂ measured to 10 bar the interval is 13 times narrower and the answer is a certain YES.', 6, fx=0.97)
    t.mark()
    t.save(26.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, default=ROOT / 'captures'); ap.add_argument('--only', default='screen,extrapolate,measure')
    args = ap.parse_args()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT / 'docs'))
    handler.log_message = lambda *a, **k: None
    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), handler); threading.Thread(target=server.serve_forever, daemon=True).start()
    base = 'http://127.0.0.1:%d' % server.server_address[1]
    with sync_playwright() as p:
        browser = p.chromium.launch(); page = browser.new_page(viewport=VIEW, device_scale_factor=1)
        if 'screen' in args.only: screen_materials(page, args.out, base)
        if 'extrapolate' in args.only: trust_extrapolation(page, args.out, base)
        if 'measure' in args.only: plan_measurement(page, args.out, base)
        browser.close()
    server.shutdown()
    print('captures written to', args.out)


if __name__ == '__main__':
    main()
