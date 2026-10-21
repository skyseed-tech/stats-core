// Shared maths for the discrete distributions (binomial, Poisson): accurate point probabilities for large n / λ,
// range sums, the exact binomial test and the normal approximation. Pure ESM, no DOM.
// Sources:
//  * Point probabilities: Loader (2000), "Fast and Accurate Computation of Binomial Probabilities" (saddle-point algorithm with the
//    Stirling error `stirlerr` and the deviance term `bd0`) — the same algorithm as R's dbinom_raw / dpois_raw (nmath), which scipy.stats.binom.pmf
//    and scipy.stats.poisson.pmf also reproduce to about 1e-14.
//  * Ranges are sums of point probabilities taken with the ratio recurrence p(k+1)/p(k) = (n−k)/(k+1)·p/(1−p) (binomial) and λ/(k+1) (Poisson),
//    starting from the largest term inside the range, so tiny tails keep full relative accuracy.
//  * Exact two-sided binomial test: sum of P(X = k) over all k with P(X = k) ≤ P(X = x)·(1 + 1e-7) — same definition and tolerance as
//    scipy.stats.binomtest(..., alternative='two-sided') and R binom.test.
//  * Normal approximation with continuity correction: Φ((b + ½ − μ)/σ) − Φ((a − ½ − μ)/σ) (e.g. Moore & McCabe, Introduction to the Practice of Statistics).
import { needNumber, needInteger, needProbability, needNonNegative } from './validate.js';
import { lnGamma, normCdf } from './stats.js';

const HALF_LN_2PI = 0.5 * Math.log(2 * Math.PI);

/** Stirling error ln n! − [ (n+½) ln n − n + ½ ln 2π ]. */
export function stirlerr(n) {
  if (n === 0) return 0;
  if (n <= 15) return lnGamma(n + 1) - (n + 0.5) * Math.log(n) + n - HALF_LN_2PI;
  const nn = n * n;
  const S0 = 1 / 12; const S1 = 1 / 360; const S2 = 1 / 1260; const S3 = 1 / 1680; const S4 = 1 / 1188;
  if (n > 500) return (S0 - S1 / nn) / n;
  if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n;
  if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

/** Deviance term bd0(x, np) = x ln(x/np) + np − x, stable when x ≈ np. */
export function bd0(x, np) {
  if (Math.abs(x - np) < 0.1 * (x + np)) {
    let v = (x - np) / (x + np);
    let s = (x - np) * v;
    if (Math.abs(s) < 1e-300) return s;
    let ej = 2 * x * v;
    v *= v;
    for (let j = 1; j < 1000; j++) {
      ej *= v;
      const s1 = s + ej / (2 * j + 1);
      if (s1 === s) return s1;
      s = s1;
    }
  }
  return x * Math.log(x / np) + np - x;
}

/** P(X = x) for X ~ Binomial(n, p), integer 0 ≤ x ≤ n. */
export function binomPmf(x, n, p) {
  needInteger(n, 'n', 0); needProbability(p, 'p'); needNumber(x, 'x');
  if (x < 0 || x > n) return 0;
  if (p <= 0) return x === 0 ? 1 : 0;
  if (p >= 1) return x === n ? 1 : 0;
  const q = 1 - p;
  if (x === 0) return Math.exp(p < 0.1 ? -bd0(n, n * q) - n * p : n * Math.log(q));
  if (x === n) return Math.exp(q < 0.1 ? -bd0(n, n * p) - n * q : n * Math.log(p));
  const lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
  const lf = Math.log(2 * Math.PI) + Math.log(x) + Math.log1p(-x / n);
  return Math.exp(lc - 0.5 * lf);
}

/** P(X = x) for X ~ Poisson(λ), integer x ≥ 0. */
export function poisPmf(x, lambda) {
  needNonNegative(lambda, 'lambda'); needNumber(x, 'x');
  if (x < 0) return 0;
  if (lambda === 0) return x === 0 ? 1 : 0;
  if (x === 0) return Math.exp(-lambda);
  return Math.exp(-stirlerr(x) - bd0(x, lambda)) / Math.sqrt(2 * Math.PI * x);
}

/** Generic range sum of a pmf with a ratio function ratio(k) = pmf(k+1)/pmf(k); `mode` is a maximiser. b may be Infinity. */
function rangeSum(a, b, pmf, ratio, mode, upperLimit) {
  a = Math.max(a, 0);
  if (b > upperLimit) b = upperLimit;
  if (a > b) return 0;
  const k0 = Math.min(Math.max(mode, a), b);
  const p0 = pmf(k0);
  if (!(p0 > 0)) return 0;
  let sum = p0;
  let term = p0;
  for (let k = k0; k < b; k++) {          // upwards
    term *= ratio(k);
    sum += term;
    if (term < 1e-17 * sum && k >= mode) break;
    if (term === 0) break;
  }
  term = p0;
  for (let k = k0; k > a; k--) {          // downwards: pmf(k−1) = pmf(k) / ratio(k−1)
    term /= ratio(k - 1);
    sum += term;
    if (term < 1e-17 * sum && k <= mode) break;
    if (term === 0) break;
  }
  return Math.min(sum, 1);
}

/** P(a ≤ X ≤ b) for Binomial(n, p). */
export function binomRange(a, b, n, p) {
  needInteger(n, 'n', 0); needProbability(p, 'p');
  if (p <= 0) return a <= 0 && b >= 0 ? 1 : 0;
  if (p >= 1) return a <= n && b >= n ? 1 : 0;
  const odds = p / (1 - p);
  return rangeSum(a, b, (k) => binomPmf(k, n, p), (k) => ((n - k) / (k + 1)) * odds, Math.floor((n + 1) * p), n);
}
/** P(X ≤ x). The smaller side is summed directly, the other is its complement. */
export function binomCdf(x, n, p) {
  needInteger(n, 'n', 0); needProbability(p, 'p'); needNumber(x, 'x');
  if (x < 0) return 0;
  if (x >= n) return 1;
  return x <= n * p ? binomRange(0, x, n, p) : 1 - binomRange(x + 1, n, n, p);
}
/** P(X > x). */
export function binomSf(x, n, p) {
  needInteger(n, 'n', 0); needProbability(p, 'p'); needNumber(x, 'x');
  if (x < 0) return 1;
  if (x >= n) return 0;
  return x <= n * p ? 1 - binomRange(0, x, n, p) : binomRange(x + 1, n, n, p);
}

/** P(a ≤ X ≤ b) for Poisson(λ); b may be Infinity. */
export function poisRange(a, b, lambda) {
  needNonNegative(lambda, 'lambda');
  if (lambda === 0) return a <= 0 && b >= 0 ? 1 : 0;
  return rangeSum(a, b, (k) => poisPmf(k, lambda), (k) => lambda / (k + 1), Math.floor(lambda), Infinity);
}
export function poisCdf(x, lambda) {
  needNonNegative(lambda, 'lambda'); needNumber(x, 'x');
  if (x < 0) return 0;
  return x <= lambda ? poisRange(0, x, lambda) : 1 - poisRange(x + 1, Infinity, lambda);
}
export function poisSf(x, lambda) {
  needNonNegative(lambda, 'lambda'); needNumber(x, 'x');
  if (x < 0) return 1;
  return x <= lambda ? 1 - poisRange(0, x, lambda) : poisRange(x + 1, Infinity, lambda);
}

/** Exact two-sided binomial test p-value (scipy.stats.binomtest / R binom.test): Σ P(X=k) over P(X=k) ≤ P(X=x)(1 + 1e-7). */
export function binomTestTwoSided(x, n, p) {
  needInteger(n, 'n', 0); needProbability(p, 'p'); needNumber(x, 'x');
  const d = binomPmf(x, n, p);
  const lim = d * (1 + 1e-7);
  let s = 0;
  for (let k = 0; k <= n; k++) {
    const pk = binomPmf(k, n, p);
    if (pk <= lim) s += pk;
  }
  return Math.min(1, s);
}

export const EVENTS = ['eq', 'le', 'lt', 'ge', 'gt', 'between'];

/** The selected event as an inclusive integer interval [lo, hi] (hi may be Infinity); lo > hi means an empty event. `max` = n for the binomial, Infinity for Poisson. */
export function eventBounds(kind, x, x2, max) {
  switch (kind) {
    case 'le': return [0, x];
    case 'lt': return [0, x - 1];
    case 'ge': return [x, max];
    case 'gt': return [x + 1, max];
    case 'between': return [x, x2];
    default: return [x, x];
  }
}

/** Normal approximation with continuity correction of P(lo ≤ X ≤ hi) for mean μ and SD σ (hi may be Infinity). */
export function normalApprox(lo, hi, mu, sigma) {
  if (!(sigma > 0) || lo > hi) return null;
  const upper = Number.isFinite(hi) ? normCdf((hi + 0.5 - mu) / sigma) : 1;
  const lower = normCdf((lo - 0.5 - mu) / sigma);
  return Math.max(0, upper - lower);
}
