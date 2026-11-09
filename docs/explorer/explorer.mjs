// Isotherm certifier: three views over real ISODB records. No remote services; all numbers come from
// data/*.json (built by scripts/build_explorer_data.py) and the lookups in numerics.mjs.
import { table, interval, decide, model, pointLnS, bestFit, answer } from './numerics.mjs';

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const SUB = { CO2: 'CO₂', N2: 'N₂', CH4: 'CH₄', C2H4: 'C₂H₄', C2H6: 'C₂H₆', H2: 'H₂', O2: 'O₂', C2H2: 'C₂H₂', C3H6: 'C₃H₆', C3H8: 'C₃H₈' };
const gas = (g) => SUB[g] || g;
const C = { ink: '#30363d', soft: '#67727d', rule: '#dfe4e8', band: '#7fb0d8', line: '#2973ac', fit: '#c0564b', fitOff: '#d9a79f', read: '#e2a03f', hinge: '#3f9f93', grey: '#9aa5ae' };
const cache = {};
const load = async (name) => (cache[name] ??= await (await fetch(`data/${name}.json`)).json());

export function fmtS(ln) {
  if (!Number.isFinite(ln)) return ln > 0 ? '∞' : '0';
  const s = Math.exp(ln);
  if (s >= 1e6) { const k = Math.floor(Math.log10(s)); return `${(s / 10 ** k).toFixed(1)}×10${String(k).replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}`; }
  return s >= 1000 ? Math.round(s).toLocaleString('en-US') : s >= 100 ? s.toFixed(0) : s >= 10 ? s.toFixed(1) : s.toFixed(2);
}
const fmtP = (p) => (p >= 100 ? p.toFixed(0) : p >= 10 ? p.toFixed(1) : p >= 1 ? p.toFixed(2) : p.toPrecision(2));
const factor = (iv) => (Number.isFinite(iv.hi - iv.lo) ? '×' + (Math.exp(iv.hi - iv.lo) >= 100 ? Math.exp(iv.hi - iv.lo).toFixed(0) : Math.exp(iv.hi - iv.lo).toFixed(1)) : 'unbounded');
const chip = (a) => `<span class="chip" data-a="${a}">${a}</span>`;
const models = (fits) => Object.fromEntries(Object.entries(fits).map(([k, v]) => [k, { ...model(k, v.theta), ok: v.ok, rmse: v.rmse }]));
const decades = (lo, hi) => { const t = []; for (let k = Math.ceil(Math.log10(lo) - 1e-9); k <= Math.floor(Math.log10(hi) + 1e-9); k++) t.push(k); return t; };
const pow10 = (k) => (k >= -2 && k <= 3 ? String(10 ** k) : `10${String(k).replace('-', '⁻').replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}`);

// ---------- plots (inline SVG)
function isoPlot(o) {
  const W = o.width || 300, H = o.height || 230, m = { l: 40, r: 10, t: 22, b: 34 };
  const x = (p) => m.l + (Math.log(Math.min(Math.max(p, o.xmin), o.xmax) / o.xmin) / Math.log(o.xmax / o.xmin)) * (W - m.l - m.r);
  const y = (q) => H - m.b - (Math.min(Math.max(q, 0), o.ymax * 1.5) / o.ymax) * (H - m.t - m.b);
  const id = 'clip' + Math.random().toString(36).slice(2, 8);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.label)}"><defs><clipPath id="${id}"><rect x="${m.l}" y="${m.t}" width="${W - m.l - m.r}" height="${H - m.t - m.b}"/></clipPath></defs>`;
  for (const k of decades(o.xmin, o.xmax)) s += `<line x1="${x(10 ** k)}" x2="${x(10 ** k)}" y1="${m.t}" y2="${H - m.b}" stroke="${C.rule}"/><text x="${x(10 ** k)}" y="${H - m.b + 14}" font-size="10" fill="${C.soft}" text-anchor="middle">${pow10(k)}</text>`;
  for (let i = 0; i <= 2; i++) { const v = (o.ymax * i) / 2; s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(v)}" y2="${y(v)}" stroke="${C.rule}"/><text x="${m.l - 5}" y="${y(v) + 3}" font-size="10" fill="${C.soft}" text-anchor="end">${v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v.toPrecision(2)}</text>`; }
  s += `<g clip-path="url(#${id})">`;
  if (o.read) s += `<rect x="${x(o.read[0])}" y="${m.t}" width="${Math.max(2, x(o.read[1]) - x(o.read[0]))}" height="${H - m.t - m.b}" fill="${C.read}" opacity=".16"/><line x1="${x(o.read[0])}" x2="${x(o.read[0])}" y1="${m.t}" y2="${H - m.b}" stroke="${C.read}" stroke-dasharray="4 3"/><line x1="${x(o.read[1])}" x2="${x(o.read[1])}" y1="${m.t}" y2="${H - m.b}" stroke="${C.read}" stroke-dasharray="4 3"/>`;
  for (const b of o.bands || []) {
    if (!b.data) continue;
    const up = b.data.p.map((p, i) => `${x(p)},${y(b.data.hi[i])}`), dn = b.data.p.map((p, i) => `${x(p)},${y(b.data.lo[i])}`).reverse();
    s += `<polygon points="${up.concat(dn).join(' ')}" fill="${b.color}" opacity="${b.opacity || 0.38}"/><polyline points="${up.join(' ')}" fill="none" stroke="${b.color}" stroke-width="1.2"/><polyline points="${dn.join(' ')}" fill="none" stroke="${b.color}" stroke-width="1.2"/>`;
  }
  for (const f of o.fits || []) {
    const pts = []; for (let i = 0; i <= 60; i++) { const p = o.xmin * (o.xmax / o.xmin) ** (i / 60); pts.push(`${x(p)},${y(f.q(p))}`); }
    s += `<polyline points="${pts.join(' ')}" fill="none" stroke="${f.ok ? C.fit : C.fitOff}" stroke-width="${f.ok ? 1.3 : 1}" opacity=".9"/>`;
  }
  if (o.vline) s += `<line x1="${x(o.vline)}" x2="${x(o.vline)}" y1="${m.t}" y2="${H - m.b}" stroke="${C.soft}" stroke-dasharray="2 3"/>`;
  for (const p of o.points || []) s += `<circle cx="${x(p.p)}" cy="${y(p.q)}" r="${p.hidden ? 2.6 : 3}" fill="${p.hidden ? '#fff' : C.ink}" stroke="${p.hidden ? (p.stroke || C.ink) : '#fff'}" stroke-width="${p.hidden ? 1.2 : 0.8}"/>`;
  s += `</g><line x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}" stroke="${C.soft}"/><text x="${m.l}" y="13" font-size="11.5" font-weight="600" fill="${C.ink}">${esc(o.title || '')}</text><text x="${(W + m.l - m.r) / 2}" y="${H - 5}" font-size="10" fill="${C.soft}" text-anchor="middle">${esc(o.xlabel || 'pressure (bar)')}</text>`;
  if (o.vline && o.vlabel) s += `<text x="${Math.min(x(o.vline) + 4, W - 70)}" y="${m.t + 10}" font-size="9.5" fill="${C.soft}">${esc(o.vlabel)}</text>`;
  return s + '</svg>';
}

function stripPlot(o) {
  // o.rows: [{label, lo, hi, color}], o.ticks: [ln S], o.best: ln S, o.tau: ln S
  const W = 620, rowH = 26, top = 30, H = top + o.rows.length * rowH + 44, l = 118, r = 20;
  const vals = [o.tau, ...(o.ticks || []).map((t) => t.v), ...(o.best != null ? [o.best] : []), ...o.rows.flatMap((d) => [d.lo, d.hi])].filter(Number.isFinite);
  let a = Math.min(...vals) - 0.5, b = Math.max(...vals) + 0.5; if (b - a < 2.4) { const c = (a + b) / 2; a = c - 1.2; b = c + 1.2; }
  const x = (v) => l + ((Math.min(Math.max(v, a), b) - a) / (b - a)) * (W - l - r);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.label)}">`;
  const L10 = Math.LN10;
  for (const k of decades(Math.exp(a), Math.exp(b))) s += `<line x1="${x(k * L10)}" x2="${x(k * L10)}" y1="${top - 6}" y2="${H - 30}" stroke="${C.rule}"/><text x="${x(k * L10)}" y="${H - 16}" font-size="10.5" fill="${C.soft}" text-anchor="middle">${pow10(k)}</text>`;
  o.rows.forEach((d, i) => {
    const yy = top + i * rowH + 9, lo = Number.isFinite(d.lo) ? x(d.lo) : l, hi = Number.isFinite(d.hi) ? x(d.hi) : W - r;
    s += `<text x="${l - 10}" y="${yy + 4}" font-size="12" fill="${C.ink}" text-anchor="end">${esc(d.label)}</text><rect x="${lo}" y="${yy - 5}" width="${Math.max(2, hi - lo)}" height="10" rx="2" fill="${d.color}"/>`;
    if (d.text) s += `<text x="${Math.min(hi + 6, W - r - 4)}" y="${yy - 8}" font-size="9.5" fill="${C.soft}" text-anchor="${hi + 90 > W ? 'end' : 'start'}">${esc(d.text)}</text>`;
  });
  const y0 = top + 9, y1 = top + (o.rows.length - 1) * rowH + 9;
  for (const t of o.ticks || []) s += `<line x1="${x(t.v)}" x2="${x(t.v)}" y1="${y1 - 9}" y2="${y1 + 9}" stroke="${t.inCap ? C.fit : C.fitOff}" stroke-width="1.6" ${t.inCap ? '' : 'stroke-dasharray="3 2"'}/>`;
  if (o.best != null && Number.isFinite(o.best)) s += `<path d="M${x(o.best)} ${y1 - 11}l5 6-5 6-5-6z" fill="${C.ink}"/>`;
  s += `<line x1="${x(o.tau)}" x2="${x(o.tau)}" y1="${y0 - 16}" y2="${H - 30}" stroke="${C.ink}" stroke-dasharray="4 3"/><text x="${x(o.tau)}" y="${y0 - 19}" font-size="11.5" fill="${C.ink}" text-anchor="middle">target S ≥ ${fmtS(o.tau)}</text>`;
  s += `<text x="${(W + l - r) / 2}" y="${H - 2}" font-size="10.5" fill="${C.soft}" text-anchor="middle">${esc(o.xlabel)}</text></svg>`;
  return s;
}

// Selectivities of all pairs of fits within tolerance; `inCap` keeps those whose Henry constants respect the caps HA, HB.
const pairEnsemble = (mA, mB, P, y, HA = Infinity, HB = Infinity) => {
  const okA = Object.values(mA).filter((m) => m.ok), okB = Object.values(mB).filter((m) => m.ok), out = [];
  for (const a of okA) for (const b of okB) { const v = pointLnS(a, b, P, y); if (Number.isFinite(v)) out.push({ v, inCap: a.henry <= HA * 1.0000001 && b.henry <= HB * 1.0000001 }); }
  return out;
};
const span = (e) => (e.length ? `${fmtS(Math.min(...e.map((d) => d.v)))}–${fmtS(Math.max(...e.map((d) => d.v)))}` : null);
const best = (side) => (side._m ??= models(side.fits))[bestFit(side.fits)];

// ---------- view 1: screen materials
const scr = { rows: [], sel: null, filter: 'ALL', q: null, data: null };
async function runScreen() {
  const f = $('query'), slug = f.mix.value, d = await load('screen-' + slug);
  const P = +f.pressure.value, y = +f.feed.value, tau = +f.tau.value, cap = f.cap.value, lnT = Math.log(tau);
  scr.data = d; scr.q = { P, y, tau, cap, slug };
  scr.rows = d.pairs.map((pr) => {
    const A = table(pr.A.r, pr.A.tables[cap]), B = table(pr.B.r, pr.B.tables[cap]);
    const iv = interval(A, B, P, y);
    return { pr, iv, ans: decide(A, B, P, y, tau), best: pointLnS(best(pr.A), best(pr.B), P, y) };
  });
  const rank = { YES: 0, UNDETERMINED: 1, NO: 2 };
  scr.rows.sort((a, b) => rank[a.ans] - rank[b.ans] || b.iv.lo - a.iv.lo);
  const n = (a) => scr.rows.filter((r) => r.ans === a).length, short = scr.rows.filter((r) => r.best >= lnT);
  $('screen-tag').textContent = `${d.pairs.length} real ${gas(d.gasA)}/${gas(d.gasB)} pairs`;
  $('screen-summary').innerHTML = `<strong>Query complete.</strong> Of ${scr.rows.length} pairs, <strong>${n('YES')}</strong> certainly meet <i>S</i> ≥ ${tau}, <strong>${n('NO')}</strong> certainly do not, and <strong>${n('UNDETERMINED')}</strong> are undetermined. A best-fit screen would shortlist ${short.length}; <strong>${short.filter((r) => r.ans !== 'YES').length}</strong> of those are not certain.`;
  $('screen-meta').textContent = `${gas(d.gasA)}/${gas(d.gasB)} · P = ${P} bar · y = ${y} · Henry cap ${cap}× minimum`;
  $('filters').innerHTML = ['ALL', 'YES', 'UNDETERMINED', 'NO'].map((a) => `<button type="button" data-f="${a}" aria-pressed="${scr.filter === a}">${a === 'ALL' ? 'All' : a.charAt(0) + a.slice(1).toLowerCase()} ${a === 'ALL' ? scr.rows.length : n(a)}</button>`).join('');
  $('export').disabled = false; $('export-status').textContent = '';
  if (!scr.sel || !scr.rows.some((r) => r.pr === scr.sel)) {     // default: a typical material that a best-fit screen would shortlist but the tables cannot decide
    const und = scr.rows.filter((r) => r.ans === 'UNDETERMINED' && r.best >= lnT && Number.isFinite(r.iv.hi - r.iv.lo)).sort((a, b) => (a.iv.hi - a.iv.lo) - (b.iv.hi - b.iv.lo));
    scr.sel = (und[Math.floor(und.length / 3)] || scr.rows[0]).pr;
  }
  drawRows(); drawPair();
}
function drawRows() {
  $('rows').innerHTML = scr.rows.filter((r) => scr.filter === 'ALL' || r.ans === scr.filter).map((r) =>
    `<tr data-id="${r.pr.id}" class="${r.pr === scr.sel ? 'selected' : ''}"><td class="name"><button type="button" title="${esc(r.pr.A.adsorbent)} · ${r.pr.A.T} K">${esc(r.pr.A.adsorbent)}</button></td><td>${fmtS(r.best)}</td><td>${fmtS(r.iv.lo)} – ${fmtS(r.iv.hi)}</td><td>${chip(r.ans)}</td></tr>`).join('');
}
function drawPair() {
  const row = scr.rows.find((r) => r.pr === scr.sel); if (!row) return;
  const { pr, iv } = row, { P, y, tau } = scr.q, d = scr.data, mA = (pr.A._m ??= models(pr.A.fits)), mB = (pr.B._m ??= models(pr.B.fits));
  $('pair-title').textContent = pr.A.adsorbent;
  $('pair-meta').innerHTML = `${pr.A.T} K · tolerance ±${(pr.A.eps * 100).toFixed(1)}% (${gas(d.gasA)}), ±${(pr.B.eps * 100).toFixed(1)}% (${gas(d.gasB)}) of full scale · source <a href="https://doi.org/${esc(pr.A.doi)}">${esc(pr.A.doi)}</a>`;
  const ens = pairEnsemble(mA, mB, P, y, pr.A.tables[scr.q.cap].H, pr.B.tables[scr.q.cap].H), inCap = ens.filter((d) => d.inCap);
  $('pair-strip').innerHTML = stripPlot({ label: 'Certified selectivity interval and model answers', tau: Math.log(tau), ticks: ens, best: row.best,
    rows: [{ label: 'Certified interval', lo: iv.lo, hi: iv.hi, color: C.line, text: factor(iv) + ' wide' }, { label: 'Fitted models', lo: NaN, hi: NaN, color: 'none' }], xlabel: `IAST selectivity ${gas(d.gasA)}/${gas(d.gasB)} at ${P} bar, y = ${y}  (◆ best fit · | fits within tolerance; dashed if beyond the Henry cap)` }).replace(/<rect[^>]*fill="none"\/>/, '');
  const panel = (side, g, read) => {
    const pm = side.p[side.p.length - 1], m = side._m, ymax = Math.max(...side.q) * 1.9;
    return isoPlot({ label: `${g} isotherm with certified band`, title: `${gas(g)} record`, xmin: side.band.p[0], xmax: side.band.p[side.band.p.length - 1], ymax,
      bands: [{ data: side.band, color: C.band }], fits: Object.values(m).filter((f) => f.ok), points: side.p.map((p, i) => ({ p, q: side.q[i] })), vline: pm, vlabel: 'last point', read });
  };
  $('pair-curves').innerHTML = `<div class="two">${panel(pr.A, d.gasA, iv.rA.slice().sort((a, b) => a - b))}${panel(pr.B, d.gasB, iv.rB.slice().sort((a, b) => a - b))}</div>`;
  const pmB = pr.B.p[pr.B.p.length - 1], far = Math.max(...iv.rB.filter(Number.isFinite), 0);
  $('screen-value').innerHTML = `<strong>Use the result.</strong> ${row.ans === 'UNDETERMINED' ? `The tables cannot decide this material: ${inCap.length > 1 ? `${inCap.length} model pairs that fit both records within tolerance and cap give selectivities of ${span(inCap)}, and` : 'the fitted models give one estimate, but'} other isotherms consistent with the same data reach ${fmtS(iv.lo)} and ${fmtS(iv.hi)}.` : `Every isotherm pair consistent with the two tables gives the same answer (${row.ans}).`} ${far > pmB * 1.05 ? `IAST reads the ${gas(d.gasB)} record at up to ${far >= 1000 ? 'thousands of' : fmtP(far)} bar, ${far / pmB >= 100 ? 'more than 100' : (far / pmB).toFixed(1)}× beyond its last measured point (${fmtP(pmB)} bar).` : `IAST reads both records inside their measured windows here.`} Export the table, or raise the Henry cap to see which answers depend on it.`;
}
function exportCSV() {
  const { P, y, tau, cap } = scr.q, d = scr.data;
  const lines = ['adsorbent,temperature_K,doi,gas_A,gas_B,P_bar,y_A,target_S,henry_cap_multiplier,best_fit_S,S_lower,S_upper,answer,record_A,record_B'];
  for (const r of scr.rows) lines.push([`"${r.pr.A.adsorbent.replace(/"/g, '""')}"`, r.pr.A.T, r.pr.A.doi, d.gasA, d.gasB, P, y, tau, cap, Math.exp(r.best).toPrecision(6), Math.exp(r.iv.lo).toPrecision(6), Number.isFinite(r.iv.hi) ? Math.exp(r.iv.hi).toPrecision(6) : 'inf', r.ans, r.pr.A.file, r.pr.B.file].join(','));
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'text/csv' }));
  a.download = `isocert_${d.gasA}_${d.gasB}_P${P}_y${y}_S${tau}.csv`; a.click(); URL.revokeObjectURL(a.href);
  $('export-status').textContent = `Saved ${scr.rows.length} rows with the query conditions.`;
}

// ---------- view 2: trust an extrapolation
const ex = { rec: null, revealed: false };
function drawExtrapolate() {
  const r = ex.rec; if (!r) return;
  const m = (r._m ??= models(r.fits)), bname = bestFit(r.fits), h = r.hidden, e = r.e, n = h.p.length;
  const okN = Object.values(r.fits).filter((f) => f.ok).length;
  const inBand = h.q.filter((q, i) => q >= h.lo[i] - e && q <= h.hi[i] + e).length, inHinge = h.q.filter((q, i) => q >= h.hlo[i] - e && q <= h.hhi[i] + e).length;
  const inBest = h.q.filter((q, i) => Math.abs(m[bname].q(h.p[i]) - q) <= e).length, last = n - 1, errLast = Math.abs(m[bname].q(h.p[last]) - h.q[last]) / h.q[last];
  $('ex-title').textContent = `${r.adsorbent} · ${gas(r.gas)} · ${r.T} K`;
  $('ex-meta').innerHTML = `Source <a href="https://doi.org/${esc(r.doi)}">${esc(r.doi)}</a>`;
  const rows = [['Visible points', `${r.p.length} (up to ${fmtP(r.p[r.p.length - 1])} bar)`], ['Tolerance', `±${e.toPrecision(2)} mmol/g (${(r.eps * 100).toFixed(1)}%)`], ['Fits within tolerance', `${okN} of ${Object.keys(r.fits).length}`],
    ['Band at the last pressure', `${h.lo[last].toPrecision(3)} – ${h.hi[last].toPrecision(3)} mmol/g`], ['Fits at the last pressure', `${Math.min(...Object.values(m).filter((f) => f.ok || okN === 0).map((f) => f.q(h.p[last]))).toPrecision(3)} – ${Math.max(...Object.values(m).filter((f) => f.ok || okN === 0).map((f) => f.q(h.p[last]))).toPrecision(3)} mmol/g`]];
  if (ex.revealed) rows.push(['Hidden points', String(n)], ['Inside the certified band', `${inBand} of ${n}`], [`Best fit (${bname}) within tolerance`, `${inBest} of ${n}`], ['Inside the shape-only band', `${inHinge} of ${n}`]);
  $('ex-facts').innerHTML = rows.map(([a, b]) => `<dt>${a}</dt><dd>${b}</dd>`).join('');
  $('ex-verdict').innerHTML = ex.revealed ? `<strong>${inBand} of ${n}</strong> hidden measurements lie inside the certified band; the best-fit model stays within tolerance on ${inBest}. At the last pressure its error is ${(errLast * 100).toFixed(0)}%.` : `${okN} model${okN === 1 ? '' : 's'} match${okN === 1 ? 'es' : ''} the visible window within tolerance. Reveal the hidden measurements to see which predictions hold.`;
  $('ex-verdict').dataset.outcome = ex.revealed && inBand < n ? 'refine' : 'ok';
  const pts = r.p.map((p, i) => ({ p, q: r.q[i] })); if (ex.revealed) h.p.forEach((p, i) => pts.push({ p, q: h.q[i], hidden: true }));
  $('ex-plot').innerHTML = isoPlot({ width: 620, height: 380, label: 'Visible points, model fits, certified band and hidden measurements', xmin: r.p[0] * 0.8, xmax: h.p[last] * 1.12, ymax: Math.max(...h.q, ...r.q) * 1.35,
    bands: [$('ex-hinge').checked ? { data: r.hinge, color: C.hinge, opacity: 0.16 } : {}, { data: r.band, color: C.band }], fits: $('ex-fits').checked ? Object.values(m) : [], points: pts, vline: r.p[r.p.length - 1], vlabel: 'visible | hidden' });
  $('ex-reveal').textContent = ex.revealed ? 'Hide measurements again' : 'Reveal hidden measurements';
}
async function initExtrapolate() {
  const d = await load('extrapolate');
  $('ex-record').innerHTML = d.records.map((r) => `<option value="${r.id}">${esc(r.adsorbent)} · ${gas(r.gas)} · ${r.T} K</option>`).join('');
  const pick = () => { ex.rec = d.records.find((r) => r.id === +$('ex-record').value); ex.revealed = false; drawExtrapolate(); };
  $('ex-record').addEventListener('change', pick); $('ex-fits').addEventListener('change', drawExtrapolate); $('ex-hinge').addEventListener('change', drawExtrapolate);
  $('ex-reveal').addEventListener('click', () => { ex.revealed = !ex.revealed; drawExtrapolate(); });
  pick();
}

// ---------- view 3: plan a measurement
const me = { pr: null };
const ivAt = (pr, i, j) => interval(table(pr.A.levels[i].r, pr.A.levels[i].table), table(pr.B.levels[j].r, pr.B.levels[j].table), pr.P, pr.y);
function drawMeasure() {
  const pr = me.pr; if (!pr) return;
  const i = +$('me-a').value, j = +$('me-b').value, tau = +$('me-tau').value, lnT = Math.log(tau), LA = pr.A.levels[i], LB = pr.B.levels[j];
  const iv = ivAt(pr, i, j), base = ivAt(pr, 0, 0), gA = gas(pr.A.gas), gB = gas(pr.B.gas);
  $('me-a-out').textContent = `${fmtP(LA.p_max)} bar · ${LA.n} pts`; $('me-b-out').textContent = `${fmtP(LB.p_max)} bar · ${LB.n} pts`;
  const mA = models(LA.fits), mB = models(LB.fits), ens = pairEnsemble(mA, mB, pr.P, pr.y, LA.table.H, LB.table.H), bst = pointLnS(mA[bestFit(LA.fits)], mB[bestFit(LB.fits)], pr.P, pr.y);
  $('me-meta').innerHTML = `${esc(pr.A.adsorbent)} · ${pr.A.T} K · ${gA}/${gB} at ${pr.P} bar, y = ${pr.y} · source <a href="https://doi.org/${esc(pr.A.doi)}">${esc(pr.A.doi)}</a>${LA.admitted && LB.admitted ? '' : ' · <strong>a record exceeds the 3% admission gate at this extent</strong>'}`;
  $('me-answer').innerHTML = `${chip(answer(iv, lnT))} &nbsp;<i>S</i> ∈ <b>[${fmtS(iv.lo)}, ${fmtS(iv.hi)}]</b> (${factor(iv)} wide).<br/>Best fit says <b>${fmtS(bst)}</b>${ens.length > 1 ? `; ${ens.length} model pairs within tolerance span <b>${span(ens)}</b>` : ''}.`;
  const seen = [], top = pr.advice.filter((a) => a.lo != null && a.hi != null).filter((a) => { const k = a.gas + Math.round(Math.log10(a.p) * 2); if (seen.includes(k)) return false; seen.push(k); return true; }).slice(0, 3);
  $('me-advice').innerHTML = top.map((a) => `<li>Measure ${a.gas === 'A' ? gA : gB} at ${fmtP(a.p)} bar<span>→ [${fmtS(a.lo)}, ${fmtS(a.hi)}]</span></li>`).join('');
  // predicted width of ln S for each candidate pressure
  const W = 470, H = 176, m = { l: 40, r: 10, t: 20, b: 32 }, adv = pr.advice.filter((a) => Number.isFinite(a.w)), wmax = Math.max(base.hi - base.lo, ...adv.map((a) => a.w)) * 1.3;
  const x = (p) => m.l + (Math.log(p / 1e-4) / Math.log(1e6)) * (W - m.l - m.r), yy = (w) => H - m.b - (w / wmax) * (H - m.t - m.b);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Predicted interval width for each candidate measurement">`;
  for (const k of [-4, -3, -2, -1, 0, 1, 2]) s += `<line x1="${x(10 ** k)}" x2="${x(10 ** k)}" y1="${m.t}" y2="${H - m.b}" stroke="${C.rule}"/><text x="${x(10 ** k)}" y="${H - m.b + 13}" font-size="10" fill="${C.soft}" text-anchor="middle">${pow10(k)}</text>`;
  for (const v of [0, wmax / 2]) s += `<text x="${m.l - 5}" y="${yy(v) + 3}" font-size="10" fill="${C.soft}" text-anchor="end">${v.toFixed(1)}</text>`;
  s += `<line x1="${m.l}" x2="${W - m.r}" y1="${yy(base.hi - base.lo)}" y2="${yy(base.hi - base.lo)}" stroke="${C.soft}" stroke-dasharray="4 3"/><text x="${W - m.r - 2}" y="${yy(base.hi - base.lo) - 6}" font-size="9.5" fill="${C.soft}" text-anchor="end">width with no new point</text>`;
  for (const [g, col] of [['A', C.fit], ['B', C.line]]) { const pts = adv.filter((a) => a.gas === g).sort((a, b) => a.p - b.p); s += `<polyline points="${pts.map((a) => `${x(a.p)},${yy(a.w)}`).join(' ')}" fill="none" stroke="${col}" stroke-width="1.5"/>` + pts.map((a) => `<circle cx="${x(a.p)}" cy="${yy(a.w)}" r="2.2" fill="${col}"/>`).join(''); }
  s += `<line x1="${x(1)}" x2="${x(1)}" y1="${m.t}" y2="${H - m.b}" stroke="${C.ink}" stroke-dasharray="2 3"/><text x="${x(1) - 4}" y="${H - m.b - 6}" font-size="9.5" fill="${C.soft}" text-anchor="end">measured ≤ 1 bar</text>`;
  s += `<line x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}" stroke="${C.soft}"/><text x="${m.l - 34}" y="${m.t - 3}" font-size="10" fill="${C.soft}">predicted width of ln S</text><text x="${(W + m.l) / 2}" y="${H - 4}" font-size="10" fill="${C.soft}" text-anchor="middle">pressure of the added point (bar) — <tspan fill="${C.fit}">${gA}</tspan> / <tspan fill="${C.line}">${gB}</tspan></text></svg>`;
  $('me-advice-plot').innerHTML = s;
  $('me-plot-meta').textContent = `Query: ${gA}/${gB} at ${pr.P} bar, y = ${pr.y}`;
  $('me-strip').innerHTML = stripPlot({ label: 'Certified interval before and after adding measurements', tau: lnT, ticks: ens, best: bst,
    rows: [{ label: 'Data up to 1 bar', lo: base.lo, hi: base.hi, color: C.grey, text: factor(base) }, { label: 'Current data', lo: iv.lo, hi: iv.hi, color: C.line, text: factor(iv) }], xlabel: `IAST selectivity ${gA}/${gB}  (◆ best fit · | fits within tolerance; dashed if beyond the Henry cap)` });
  const panel = (side, L, g, read) => isoPlot({ label: `${g} record`, title: `${gas(g)} record`, xmin: L.band.p[0], xmax: L.band.p[L.band.p.length - 1], ymax: Math.max(...side.q) * 1.35, bands: [{ data: L.band, color: C.band }],
    points: side.p.map((p, k) => ({ p, q: side.q[k], hidden: k >= L.n, stroke: C.grey })), vline: L.p_max, vlabel: 'last point in use', read });
  $('me-curves').innerHTML = `<div class="two">${panel(pr.A, LA, pr.A.gas, iv.rA.slice().sort((a, b) => a - b))}${panel(pr.B, LB, pr.B.gas, iv.rB.slice().sort((a, b) => a - b))}</div>`;
  const nA = pr.A.levels.length - 1, nB = pr.B.levels.length - 1, onlyA = ivAt(pr, nA, 0), onlyB = ivAt(pr, 0, nB);
  $('me-value').innerHTML = i + j === 0
    ? `<strong>Start here.</strong> With both records cut at 1 bar the interval is ${factor(base)} wide and the answer to <i>S</i> ≥ ${tau} is ${answer(base, lnT)}. The advice ranks a ${top[0] ? (top[0].gas === 'A' ? gA : gB) + ' point near ' + fmtP(top[0].p) + ' bar' : 'new point'} first. Move a slider to add the real measurements.`
    : `<strong>What the measurements bought.</strong> Interval width went from ${factor(base)} to ${factor(iv)}; the answer to <i>S</i> ≥ ${tau} is now ${answer(iv, lnT)}. Extending only ${gA} to ${fmtP(pr.A.levels[nA].p_max)} bar gives ${factor(onlyA)}; extending only ${gB} to ${fmtP(pr.B.levels[nB].p_max)} bar gives ${factor(onlyB)}.`;
}
async function initMeasure() {
  const d = await load('measure');
  $('me-pair').innerHTML = d.pairs.map((p) => `<option value="${p.id}">${esc(p.A.adsorbent)} · ${gas(p.A.gas)}/${gas(p.B.gas)} · ${p.A.T} K</option>`).join('');
  const pick = () => { const pr = (me.pr = d.pairs.find((p) => p.id === +$('me-pair').value)); $('me-la').textContent = gas(pr.A.gas); $('me-lb').textContent = gas(pr.B.gas); $('me-a').max = pr.A.levels.length - 1; $('me-b').max = pr.B.levels.length - 1; $('me-a').value = 0; $('me-b').value = 0; const b0 = ivAt(pr, 0, 0); $('me-tau').value = Number.isFinite(b0.lo + b0.hi) ? +Math.exp((b0.lo + b0.hi) / 2).toPrecision(2) : pr.tau; drawMeasure(); };
  $('me-pair').addEventListener('change', pick); for (const id of ['me-a', 'me-b', 'me-tau']) $(id).addEventListener('input', drawMeasure);
  pick();
}

// ---------- shell
function route() {
  const view = ['screen', 'extrapolate', 'measure'].includes(location.hash.slice(1)) ? location.hash.slice(1) : 'screen';
  for (const s of document.querySelectorAll('main > section')) s.hidden = s.id !== view;
  for (const a of document.querySelectorAll('nav a')) a.getAttribute('href') === '#' + view ? a.setAttribute('aria-current', 'page') : a.removeAttribute('aria-current');
}
async function main() {
  const f = $('query');
  f.addEventListener('submit', (e) => { e.preventDefault(); runScreen(); });
  f.mix.addEventListener('change', async () => { const d = await load('screen-' + f.mix.value); f.feed.value = d.y; f.tau.value = d.tau; $('feed-label').textContent = gas(d.gasA); scr.sel = null; scr.filter = 'ALL'; runScreen(); });
  $('filters').addEventListener('click', (e) => { const b = e.target.closest('button'); if (!b) return; scr.filter = b.dataset.f; for (const x of $('filters').children) x.setAttribute('aria-pressed', x === b); drawRows(); });
  $('rows').addEventListener('click', (e) => { const tr = e.target.closest('tr'); if (!tr) return; scr.sel = scr.rows.find((r) => r.pr.id === +tr.dataset.id).pr; drawRows(); drawPair(); });
  $('export').addEventListener('click', exportCSV);
  const done = await Promise.allSettled([runScreen(), initExtrapolate(), initMeasure()]), failed = done.filter((d) => d.status === 'rejected');
  if (failed.length) $('load').textContent = 'Some explorer data could not be loaded: ' + failed[0].reason.message; else $('load').hidden = true;
  window.addEventListener('hashchange', route); route();
}
if (typeof document !== 'undefined') main().catch((err) => { $('load').textContent = 'Could not load the explorer data: ' + err.message; });
