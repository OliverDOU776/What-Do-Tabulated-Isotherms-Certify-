// The browser engine must reproduce the Python engine on the stored tables and fits.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { table, interval, decide, model, pointLnS, bestFit } from '../docs/explorer/numerics.mjs';

const data = (name) => JSON.parse(readFileSync(new URL(`../docs/explorer/data/${name}.json`, import.meta.url)));
const ref = data('reference');
const num = (v) => (v === 'inf' ? Infinity : v === '-inf' ? -Infinity : v);
const close = (a, b, tol) => (a === b) || Math.abs(a - b) <= tol;

test('screening intervals, answers and best-fit selectivities agree with Python', () => {
  const files = {};
  let n = 0;
  for (const c of ref.cases.filter((x) => x.view === 'screen')) {
    const d = (files[c.file] ??= data('screen-' + c.file)), pr = d.pairs.find((p) => p.id === c.id);
    const A = table(pr.A.r, pr.A.tables[c.cap]), B = table(pr.B.r, pr.B.tables[c.cap]), iv = interval(A, B, c.P, c.y);
    assert.ok(close(iv.lo, num(c.lo), 1e-9) && close(iv.hi, num(c.hi), 1e-9), `interval of pair ${c.id} in ${c.file} at P=${c.P}`);
    assert.equal(decide(A, B, c.P, c.y, c.tau), c.answer);
    const best = pointLnS(model(bestFit(pr.A.fits), pr.A.fits[bestFit(pr.A.fits)].theta), model(bestFit(pr.B.fits), pr.B.fits[bestFit(pr.B.fits)].theta), c.P, c.y);
    assert.ok(close(best, c.best, 2e-4), `best-fit ln S of pair ${c.id} in ${c.file}: ${best} vs ${c.best}`);
    n++;
  }
  assert.ok(n >= 60);
});

test('measurement-planning intervals agree with Python', () => {
  const d = data('measure');
  for (const c of ref.cases.filter((x) => x.view === 'measure')) {
    const pr = d.pairs.find((p) => p.id === c.id), LA = pr.A.levels[c.i], LB = pr.B.levels[c.j];
    const iv = interval(table(LA.r, LA.table), table(LB.r, LB.table), pr.P, pr.y);
    assert.ok(close(iv.lo, c.lo, 1e-9) && close(iv.hi, c.hi, 1e-9));
  }
});

test('threshold predicate is consistent with the interval', () => {
  const d = data('measure'), pr = d.pairs[0], LA = pr.A.levels[0], LB = pr.B.levels[0];
  const A = table(LA.r, LA.table), B = table(LB.r, LB.table), iv = interval(A, B, pr.P, pr.y);
  assert.equal(decide(A, B, pr.P, pr.y, Math.exp(iv.lo - 0.01)), 'YES');
  assert.equal(decide(A, B, pr.P, pr.y, Math.exp(iv.hi + 0.01)), 'NO');
  assert.equal(decide(A, B, pr.P, pr.y, Math.exp((iv.lo + iv.hi) / 2)), 'UNDETERMINED');
});
