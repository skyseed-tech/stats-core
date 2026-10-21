// Statistics helpers: normal, chi-square, t and F distributions, quantiles, number parsing. Pure ESM, no DOM, deterministic.
// No strings here: errors are thrown by the callers; functions return NaN / null for undefined results.

/** Arithmetic mean. */
export const mean = (xs) => xs.reduce((a, b) => a + b, 0) / xs.length;

/** Sample standard deviation (n − 1 in the denominator), computed in two passes. NaN for fewer than 2 values. */
export function sampleSD(xs) {
  const n = xs.length;
  if (n < 2) return NaN;
  const m = mean(xs);
  return Math.sqrt(xs.reduce((a, x) => a + (x - m) ** 2, 0) / (n - 1));
}

/** Population standard deviation (n in the denominator). */
export function popSD(xs) {
  const n = xs.length;
  if (n < 1) return NaN;
  const m = mean(xs);
  return Math.sqrt(xs.reduce((a, x) => a + (x - m) ** 2, 0) / n);
}

// ---- standard normal distribution -------------------------------------------------------------------------------
// erfc(z) for z >= 0:
//  * z < 2: power series erf(z) = 2/√π · e^(−z²) · Σ 2^n z^(2n+1) / (1·3·5···(2n+1))   (Abramowitz & Stegun 7.1.6),
//    erfc = 1 − erf; all terms are positive, so there is no cancellation inside the sum.
//  * z >= 2: continued fraction erfc(z) = e^(−z²)/√π · 1/(z + (1/2)/(z + 1/(z + (3/2)/(z + 2/(z + …))))) (A&S 7.1.14),
//    evaluated from the bottom up; keeps full relative accuracy in the far tail (needed for ppm of 6-sigma processes).
// Accuracy about 1e-15 (absolute) and 1e-13 (relative in the tails), far better than the ≥ 1e-7 required.
function erfcPos(z) {
  if (z < 2) {
    let term = z;
    let sum = z;
    for (let n = 1; n < 200; n++) {
      term *= (2 * z * z) / (2 * n + 1);
      sum += term;
      if (term < sum * 1e-17) break;
    }
    return 1 - (2 / Math.sqrt(Math.PI)) * Math.exp(-z * z) * sum;
  }
  let f = z;
  for (let k = 80; k >= 1; k--) f = z + (k / 2) / f;
  return Math.exp(-z * z) / Math.sqrt(Math.PI) / f;
}

/** Standard normal cumulative distribution Φ(x). */
export function normCdf(x) {
  if (Number.isNaN(x)) return NaN;
  if (x === Infinity) return 1;
  if (x === -Infinity) return 0;
  const z = Math.abs(x) / Math.SQRT2;
  const tail = 0.5 * erfcPos(z); // P(Z > |x|)
  return x >= 0 ? 1 - tail : tail;
}

/** Upper tail 1 − Φ(x) with full relative accuracy for large x (use this for ppm / DPMO). */
export const normSf = (x) => normCdf(-x);

/** Standard normal density φ(x). */
export const normPdf = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);

// Inverse Φ: Peter J. Acklam's rational approximation (relative error < 1.15e-9, https://web.archive.org/web/2015*/home.online.no/~pjacklam/notes/invnorm/),
// followed by one step of Halley refinement against normCdf, which brings it to about machine precision.
const A = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
const B = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
const C = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
const D = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];

/** Inverse standard normal Φ⁻¹(p) for 0 < p < 1; ±Infinity at 0 / 1; NaN outside [0, 1]. */
export function normInv(p) {
  if (!(p >= 0 && p <= 1)) return NaN;
  if (p === 0) return -Infinity;
  if (p === 1) return Infinity;
  const plow = 0.02425;
  let x;
  if (p < plow) {
    const q = Math.sqrt(-2 * Math.log(p));
    x = (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1);
  } else if (p > 1 - plow) {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    x = -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1);
  } else {
    const q = p - 0.5;
    const r = q * q;
    x = ((((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q) / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1);
  }
  // one Halley step (upper-tail formulation for p > 0.5 keeps relative accuracy)
  const e = p > 0.5 ? -(normSf(x) - (1 - p)) : normCdf(x) - p;
  const u = e * Math.sqrt(2 * Math.PI) * Math.exp((x * x) / 2);
  return x - u / (1 + (x * u) / 2);
}

/** Inverse of the upper tail: x such that 1 − Φ(x) = q (q = defect probability). Accurate for tiny q such as 3.4e-6. */
export const normInvSf = (q) => -normInv(q);

// ---- F distribution (for the Gage R&R interaction test) -----------------------------------------------------------
// ln Γ(x): Lanczos approximation (g = 7, n = 9), Numerical Recipes / Godfrey coefficients.
const LG = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
export function lnGamma(x) {
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lnGamma(1 - x);
  x -= 1;
  let a = LG[0];
  const tt = x + 7.5;
  for (let i = 1; i < 9; i++) a += LG[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(tt) - tt + Math.log(a);
}

// Regularized incomplete beta I_x(a, b): continued fraction of Numerical Recipes "betacf" (modified Lentz).
function betacf(a, b, x) {
  const TINY = 1e-300;
  let c = 1;
  let d = 1 - ((a + b) * x) / (a + 1);
  if (Math.abs(d) < TINY) d = TINY;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 300; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((a + m2 - 1) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < TINY) d = TINY;
    c = 1 + aa / c; if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d; h *= d * c;
    aa = (-(a + m) * (a + b + m) * x) / ((a + m2) * (a + m2 + 1));
    d = 1 + aa * d; if (Math.abs(d) < TINY) d = TINY;
    c = 1 + aa / c; if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-15) break;
  }
  return h;
}
export function betaInc(a, b, x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? (bt * betacf(a, b, x)) / a : 1 - (bt * betacf(b, a, 1 - x)) / b;
}

/** Upper-tail probability P(F > f) of the F distribution with d1, d2 degrees of freedom (the ANOVA p-value). */
export function fSf(f, d1, d2) {
  if (!(f >= 0) || !(d1 > 0) || !(d2 > 0)) return NaN;
  if (f === 0) return 1;
  if (f === Infinity) return 0;
  return betaInc(d2 / 2, d1 / 2, d2 / (d2 + d1 * f));
}

// ---- SPC constants ---------------------------------------------------------------------------------------------
// d2: expected value of the range of n standard-normal observations. Published SPC constants (AIAG SPC manual 2nd ed.;
// Montgomery, Introduction to Statistical Quality Control, Appendix VI). σ̂ = R̄ / d2.
export const D2 = { 2: 1.128, 3: 1.693, 4: 2.059, 5: 2.326, 6: 2.534, 7: 2.704, 8: 2.847, 9: 2.970, 10: 3.078 };

// ---- parsing -----------------------------------------------------------------------------------------------------
/**
 * Numbers from free text: separated by commas, semicolons, spaces, tabs or new lines. A decimal comma is not supported.
 * Returns { values: number[], bad: string[] } — `bad` lists tokens that are not numbers (the caller decides what to do).
 */
export function parseNumberList(text) {
  const values = [];
  const bad = [];
  for (const raw of String(text ?? '').split(/[\s,;]+/)) {
    if (!raw) continue;
    const s = raw.replace(/−/g, '-');
    if (/^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/i.test(s)) values.push(Number(s));
    else bad.push(raw);
  }
  return { values, bad };
}

/** One number from a typed field ("" → NaN, "1,000" → 1000). */
export function parseOne(v) {
  if (typeof v === 'number') return v;
  const s = String(v ?? '').trim().replace(/,/g, '').replace(/−/g, '-');
  return /^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/i.test(s) ? Number(s) : NaN;
}

// ====================================================================================================================
// Chi-square, t, F distributions and quantile helpers.
// ====================================================================================================================

// ---- regularized incomplete gamma (chi-square distribution) ---------------------------------------------------------
// P(a, x) = γ(a, x)/Γ(a): series for x < a + 1, continued fraction (modified Lentz) for the complement otherwise.
// Numerical Recipes 3rd ed. §6.2 ("gser" / "gcf"); accuracy about 1e-14.
function gammaSeries(a, x) {
  let ap = a;
  let del = 1 / a;
  let sum = del;
  for (let n = 0; n < 1000; n++) {
    ap += 1;
    del *= x / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * 1e-16) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - lnGamma(a));
}
function gammaCF(a, x) {
  const TINY = 1e-300;
  let b = x + 1 - a;
  let c = 1 / TINY;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < TINY) d = TINY;
    c = b + an / c; if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return Math.exp(-x + a * Math.log(x) - lnGamma(a)) * h;
}
/** Regularized lower incomplete gamma P(a, x). NaN for a <= 0 or x < 0. */
export function gammaP(a, x) {
  if (!(a > 0) || !(x >= 0)) return NaN;
  if (x === 0) return 0;
  if (x === Infinity) return 1;
  return x < a + 1 ? gammaSeries(a, x) : 1 - gammaCF(a, x);
}
/** Regularized upper incomplete gamma Q(a, x) = 1 − P(a, x), accurate in the far tail. */
export function gammaQ(a, x) {
  if (!(a > 0) || !(x >= 0)) return NaN;
  if (x === 0) return 1;
  if (x === Infinity) return 0;
  return x < a + 1 ? 1 - gammaSeries(a, x) : gammaCF(a, x);
}

// ---- generic monotone inverse (bisection) ---------------------------------------------------------------------------
// Finds x >= 0 with fn(x) = target for a monotone fn on [0, ∞); `increasing` tells the direction. About 1e-15 relative.
function invMono(fn, target, increasing, hi0) {
  const reached = (x) => (increasing ? fn(x) >= target : fn(x) <= target);
  let lo = 0;
  let hi = hi0;
  for (let i = 0; i < 2000 && !reached(hi); i++) { lo = hi; hi *= 2; }
  for (let i = 0; i < 400; i++) {
    const mid = (lo + hi) / 2;
    if (reached(mid)) hi = mid; else lo = mid;
    if (hi - lo <= 1e-15 * hi) break;
  }
  return (lo + hi) / 2;
}

// ---- chi-square distribution (k degrees of freedom) -----------------------------------------------------------------
/** Chi-square cumulative distribution P(X <= x). */
export function chiSqCdf(x, k) {
  if (!(k > 0) || Number.isNaN(x)) return NaN;
  return x <= 0 ? 0 : gammaP(k / 2, x / 2);
}
/** Chi-square upper tail P(X > x) (the usual p-value). */
export function chiSqSf(x, k) {
  if (!(k > 0) || Number.isNaN(x)) return NaN;
  return x <= 0 ? 1 : gammaQ(k / 2, x / 2);
}
/** Chi-square density. */
export function chiSqPdf(x, k) {
  if (!(k > 0) || x < 0 || Number.isNaN(x)) return x < 0 ? 0 : NaN;
  if (x === 0) return k < 2 ? Infinity : k === 2 ? 0.5 : 0;
  return Math.exp((k / 2 - 1) * Math.log(x) - x / 2 - (k / 2) * Math.LN2 - lnGamma(k / 2));
}
/** Chi-square quantile: x with P(X <= x) = p (lower-tail critical value). */
export function chiSqInv(p, k) {
  if (!(p >= 0 && p <= 1) || !(k > 0)) return NaN;
  if (p === 0) return 0;
  if (p === 1) return Infinity;
  return invMono((x) => chiSqCdf(x, k), p, true, k + 10);
}
/** Chi-square upper critical value: x with P(X > x) = q (for the right-tail test at level q). */
export function chiSqInvSf(q, k) {
  if (!(q >= 0 && q <= 1) || !(k > 0)) return NaN;
  if (q === 1) return 0;
  if (q === 0) return Infinity;
  return invMono((x) => chiSqSf(x, k), q, false, k + 10);
}

// ---- Student t distribution (df may be fractional, e.g. Welch) --------------------------------------------------------
/** Student t upper tail P(T > t). */
export function tSf(t, df) {
  if (!(df > 0) || Number.isNaN(t)) return NaN;
  if (t === Infinity) return 0;
  if (t === -Infinity) return 1;
  const half = 0.5 * betaInc(df / 2, 0.5, df / (df + t * t));
  return t >= 0 ? half : 1 - half;
}
/** Student t cumulative distribution P(T <= t). */
export const tCdf = (t, df) => tSf(-t, df);
/** Student t density. */
export function tPdf(t, df) {
  if (!(df > 0) || Number.isNaN(t)) return NaN;
  return Math.exp(lnGamma((df + 1) / 2) - lnGamma(df / 2) - 0.5 * Math.log(df * Math.PI) - ((df + 1) / 2) * Math.log(1 + (t * t) / df));
}
/** Student t quantile: t with P(T <= t) = p, e.g. tInv(0.975, 10) = 2.228 (two-sided 95 % critical value). */
export function tInv(p, df) {
  if (!(p >= 0 && p <= 1) || !(df > 0)) return NaN;
  if (p === 0) return -Infinity;
  if (p === 1) return Infinity;
  if (p === 0.5) return 0;
  const pos = (q) => invMono((x) => tSf(x, df), q, false, 2); // x >= 0 with upper tail q
  return p < 0.5 ? -pos(p) : pos(1 - p);
}
/** Upper-tail critical value: t with P(T > t) = q. */
export const tInvSf = (q, df) => -tInv(q, df);

// ---- F distribution -------------------------------------------------------------------------------------------------
/** F cumulative distribution P(F <= f). */
export function fCdf(f, d1, d2) {
  if (!(f >= 0) || !(d1 > 0) || !(d2 > 0)) return Number.isNaN(f) || f < 0 ? (f < 0 ? 0 : NaN) : NaN;
  if (f === Infinity) return 1;
  return betaInc(d1 / 2, d2 / 2, (d1 * f) / (d1 * f + d2));
}
/** F density. */
export function fPdf(f, d1, d2) {
  if (!(d1 > 0) || !(d2 > 0) || Number.isNaN(f)) return NaN;
  if (f < 0) return 0;
  if (f === 0) return d1 < 2 ? Infinity : d1 === 2 ? 1 : 0;
  const lb = lnGamma(d1 / 2) + lnGamma(d2 / 2) - lnGamma((d1 + d2) / 2);
  return Math.exp(0.5 * (d1 * Math.log(d1 * f) + d2 * Math.log(d2) - (d1 + d2) * Math.log(d1 * f + d2)) - Math.log(f) - lb);
}
/** F upper critical value: f with P(F > f) = q (e.g. q = 0.05 gives the usual 5 % table). */
export function fInvSf(q, d1, d2) {
  if (!(q >= 0 && q <= 1) || !(d1 > 0) || !(d2 > 0)) return NaN;
  if (q === 1) return 0;
  if (q === 0) return Infinity;
  return invMono((x) => fSf(x, d1, d2), q, false, 2);
}
/** F quantile: f with P(F <= f) = p. */
export function fInv(p, d1, d2) {
  if (!(p >= 0 && p <= 1) || !(d1 > 0) || !(d2 > 0)) return NaN;
  if (p === 0) return 0;
  if (p === 1) return Infinity;
  return invMono((x) => fCdf(x, d1, d2), p, true, 2);
}

// ---- quartiles and quantiles ------------------------------------------------------------------------------------------
// There is no single universal definition (NIST/SEMATECH e-Handbook 7.2.6.2; Hyndman & Fan 1996, American Statistician 50:361).
//   halves  – median of the lower / upper half of the sorted data, the middle value left out when n is odd
//             (the textbook method of Moore & McCabe and most graphing calculators)
//   tukey   – Tukey's hinges: as 'halves' but the middle value belongs to both halves when n is odd (Tukey 1977, EDA)
//   nist    – interpolate at position p·(n + 1)   (NIST Handbook default = Hyndman & Fan type 6; "exclusive" in spreadsheets)
//   linear  – interpolate at position 1 + p·(n − 1) (Hyndman & Fan type 7, the default of R and "inclusive" in spreadsheets)
export const QUARTILE_METHODS = ['halves', 'tukey', 'nist', 'linear'];

/** Median of an already sorted array (NaN when empty). */
export function medianSorted(s) {
  const n = s.length;
  if (!n) return NaN;
  return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
}

/** Quantile of a sorted array with the 'nist' (type 6) or 'linear' (type 7) rule. */
export function quantileSorted(s, p, method = 'linear') {
  const n = s.length;
  if (!n) return NaN;
  if (method === 'nist') {
    const pos = p * (n + 1);
    const k = Math.floor(pos);
    if (k < 1) return s[0];
    if (k >= n) return s[n - 1];
    return s[k - 1] + (pos - k) * (s[k] - s[k - 1]);
  }
  const h = (n - 1) * p;
  const lo = Math.floor(h);
  return lo >= n - 1 ? s[n - 1] : s[lo] + (h - lo) * (s[lo + 1] - s[lo]);
}

/** { q1, median, q3 } of a sorted array by one of QUARTILE_METHODS. */
export function quartilesSorted(s, method = 'halves') {
  const n = s.length;
  if (!n) return { q1: NaN, median: NaN, q3: NaN };
  const median = medianSorted(s);
  if (method === 'nist' || method === 'linear') return { q1: quantileSorted(s, 0.25, method), median, q3: quantileSorted(s, 0.75, method) };
  if (n === 1) return { q1: s[0], median, q3: s[0] };
  const half = method === 'tukey' ? Math.ceil(n / 2) : Math.floor(n / 2);
  return { q1: medianSorted(s.slice(0, half)), median, q3: medianSorted(s.slice(n - half)) };
}
