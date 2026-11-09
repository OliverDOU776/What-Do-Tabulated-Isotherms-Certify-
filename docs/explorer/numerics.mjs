// Online stage of IsoCert in the browser: lookups in compiled envelope tables and binary IAST on them.
// Mirrors src/isocert/compiled.py and src/isocert/iast.py; tests/test_explorer.mjs checks agreement with Python.

export function table(r, t) {
  // r: horizon nodes; t: {U, L, q (= q_lo at the nodes), H (Henry cap), qt (= q_hi at the last node)}
  const lnr = r.map(Math.log), m = r.length - 1;
  const index = (x) => { let a = 0, b = m; while (b - a > 1) { const c = (a + b) >> 1; if (r[c] <= x) a = c; else b = c; } return a; };
  return {
    U(x) {
      if (x <= r[0]) return Math.min(t.H * x, t.U[0]);
      if (x >= r[m]) return t.U[m] + t.qt * (x / r[m] - 1);
      const j = index(x), w = (Math.log(x) - lnr[j]) / (lnr[j + 1] - lnr[j]);
      return (1 - w) * t.U[j] + w * t.U[j + 1];
    },
    L(x) {
      if (x <= r[0]) return t.L[0] * x / r[0];
      const j = x >= r[m] ? m : index(x);
      return t.L[j] + t.q[j] * (Math.log(x) - lnr[j]);
    },
  };
}

// Root of Pi1(P y (1+e^-z)) - Pi2(P (1-y) (1+e^z)), decreasing in z = ln(x/(1-x)).
export function solveZ(Pi1, Pi2, P, y, tol = 1e-3, zlo = -30, zhi = 30, side = null) {
  const f = (z) => Pi1(P * y * (1 + Math.exp(-z))) - Pi2(P * (1 - y) * (1 + Math.exp(z)));
  if (f(zlo) <= 0) return side === 'lo' ? -Infinity : zlo;
  if (f(zhi) >= 0) return side === 'up' ? Infinity : zhi;
  while (zhi - zlo > tol) { const mid = 0.5 * (zlo + zhi); if (f(mid) > 0) zlo = mid; else zhi = mid; }
  return side === 'up' ? zhi : side === 'lo' ? zlo : 0.5 * (zlo + zhi);
}

// Certified interval of ln S for envelope objects A, B (each with U and L).
export function interval(A, B, P, y) {
  const c = Math.log(y / (1 - y));
  const zu = solveZ(A.U, B.L, P, y, 1e-3, -30, 30, 'up'), zl = solveZ(A.L, B.U, P, y, 1e-3, -30, 30, 'lo');
  const rA = (z) => P * y * (1 + Math.exp(-z)), rB = (z) => P * (1 - y) * (1 + Math.exp(z));
  return { lo: zl - c, hi: zu - c, rA: [rA(zu), rA(zl)], rB: [rB(zl), rB(zu)] };
}

// Threshold predicate: four envelope values, no root solve.
export function decide(A, B, P, y, tau) {
  const z = Math.log(tau) + Math.log(y / (1 - y));
  const ra = P * y * (1 + Math.exp(-z)), rb = P * (1 - y) * (1 + Math.exp(z));
  if (A.L(ra) >= B.U(rb)) return 'YES';
  if (A.U(ra) < B.L(rb)) return 'NO';
  return 'UNDETERMINED';
}

// ---- parametric families (same parametrisation as src/isocert/baselines.py)
const sigmoid = (v) => 1 / (1 + Math.exp(-v));
function numericPi(q, r) {            // int_0^r q(t)/t dt: analytic prefix (q ~ C t near 0) + Simpson in u = ln t
  const t0 = r * 1e-12, n = 240, a = Math.log(t0), h = (Math.log(r) - a) / n;
  let s = q(t0) + q(r);
  for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * q(Math.exp(a + i * h));
  return q(t0) + s * h / 3;
}
// Each model exposes q(p), Pi(r) and its Henry constant (Infinity for Freundlich and for Sips with exponent below one).
export function model(name, th) {
  const e = th.map(Math.exp);
  if (name === 'Langmuir') return { henry: e[0] * e[1], q: (p) => e[0] * e[1] * p / (1 + e[1] * p), Pi: (r) => e[0] * Math.log1p(e[1] * r) };
  if (name === 'DSL') return { henry: e[0] * e[1] + e[2] * e[3], q: (p) => e[0] * e[1] * p / (1 + e[1] * p) + e[2] * e[3] * p / (1 + e[3] * p), Pi: (r) => e[0] * Math.log1p(e[1] * r) + e[2] * Math.log1p(e[3] * r) };
  if (name === 'Sips') { const n = sigmoid(th[2]); return { henry: n < 1 ? Infinity : e[0] * e[1], q: (p) => { const x = (e[1] * p) ** n; return e[0] * x / (1 + x); }, Pi: (r) => e[0] / n * Math.log1p((e[1] * r) ** n) }; }
  if (name === 'Toth') { const t = sigmoid(th[2]); const q = (p) => e[0] * e[1] * p / (1 + (e[1] * p) ** t) ** (1 / t); return { henry: e[0] * e[1], q, Pi: (r) => numericPi(q, r) }; }
  if (name === 'UNILAN') { const s = Math.min(Math.max(e[2], 1e-3), 30); const q = (p) => e[0] / (2 * s) * (Math.log1p(e[1] * Math.exp(s) * p) - Math.log1p(e[1] * Math.exp(-s) * p)); return { henry: e[0] * e[1] * Math.sinh(s) / s, q, Pi: (r) => numericPi(q, r) }; }
  if (name === 'Freundlich') { const a = sigmoid(th[1]); return { henry: Infinity, q: (p) => e[0] * p ** a, Pi: (r) => e[0] * r ** a / a }; }
  throw new Error('unknown family ' + name);
}

// ln S of two single isotherms.
export function pointLnS(mA, mB, P, y) {
  return solveZ(mA.Pi, mB.Pi, P, y, 1e-6, -60, 60) - Math.log(y / (1 - y));
}

export const bestFit = (fits) => Object.keys(fits).reduce((a, b) => (fits[a].rmse <= fits[b].rmse ? a : b));
export const answer = (iv, lnTau) => (iv.lo >= lnTau ? 'YES' : iv.hi < lnTau ? 'NO' : 'UNDETERMINED');
