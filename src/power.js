// Power and sample-size maths for t tests, two proportions and correlation. Pure ESM, no DOM.
// Sources:
//  * Noncentral t cumulative distribution: Lenth, R.V. (1989) "Algorithm AS 243: Cumulative distribution function of the non-central t
//    distribution", Applied Statistics 38(1), 185-189 (the same algorithm as R pt(ncp=) and scipy.stats.nct.cdf for moderate ncp).
//    For |ncp| > 37.6 or df > 4e5 the normal approximation of Abramowitz & Stegun 26.7.10 is used (as R does).
//  * t-test power: Cohen (1988) Statistical Power Analysis for the Behavioral Sciences, ch. 2 and 3 — power = P(|T'| > t_crit),
//    T' noncentral t with ncp = d·sqrt(n/2) (two independent groups of n) or d·sqrt(n) (one sample / paired). Same as statsmodels
//    TTestIndPower / TTestPower .power() and G*Power-class exact methods.
//  * Two proportions: Cohen's h = 2·asin(sqrt(p1)) − 2·asin(sqrt(p2)) (Cohen 1988 ch. 6), normal approximation,
//    ncp = h·sqrt(n/2) per group — statsmodels NormalIndPower with proportion_effectsize.
//  * Correlation: Fisher z transform, z = atanh(r), standard error 1/sqrt(n − 3) (Cohen 1988 ch. 3; Fisher 1921).
import { needNumber, needInteger, needOpenProbability } from './validate.js';
import { normCdf, normInv, tCdf, tInv, betaInc, lnGamma } from './stats.js';

const SQRT_2_OVER_PI = Math.sqrt(2 / Math.PI);
const LN_SQRT_PI = 0.5 * Math.log(Math.PI);

/** Noncentral t CDF P(T' ≤ x) with df degrees of freedom and noncentrality ncp (AS 243). */
export function nctCdf(x, df, ncp) {
  if (!(df > 0) || Number.isNaN(x) || Number.isNaN(ncp)) return NaN;
  if (ncp === 0) return tCdf(x, df);
  if (x === Infinity) return 1;
  if (x === -Infinity) return 0;
  let negdel; let tt; let del;
  if (x >= 0) { negdel = false; tt = x; del = ncp; } else {
    if (ncp > 40 && (df > 1e5 || ncp * ncp > 2 * Math.LN2 * 1021)) return 0;
    negdel = true; tt = -x; del = -ncp;
  }
  if (df > 4e5 || del * del > 2 * Math.LN2 * 1021) {
    const s = 1 / (4 * df);
    const lower = normCdf((tt * (1 - s) - del) / Math.sqrt(1 + tt * tt * 2 * s));
    return negdel ? 1 - lower : lower;
  }
  const xx0 = tt * tt;
  const xx = xx0 / (xx0 + df);
  let tnc;
  if (xx > 0) {
    const lambda = del * del;
    let p = 0.5 * Math.exp(-0.5 * lambda);
    let q = SQRT_2_OVER_PI * p * del;
    let s = 0.5 - p;
    if (s < 1e-7) s = -0.5 * Math.expm1(-0.5 * lambda);
    let a = 0.5;
    const b = 0.5 * df;
    const rxb = Math.pow(1 - xx, b);
    const albeta = LN_SQRT_PI + lnGamma(b) - lnGamma(0.5 + b);
    let xodd = betaInc(0.5, b, xx);
    let godd = 2 * rxb * Math.exp(a * Math.log(xx) - albeta);
    tnc = b * xx;
    let xeven = tnc < Number.EPSILON ? tnc : 1 - rxb;
    let geven = tnc * rxb;
    tnc = p * xodd + q * xeven;
    for (let it = 1; it <= 20000; it++) {
      a += 1;
      xodd -= godd;
      xeven -= geven;
      godd *= xx * (a + b - 1) / a;
      geven *= xx * (a + b - 0.5) / (a + 0.5);
      p *= lambda / (2 * it);
      q *= lambda / (2 * it + 1);
      tnc += p * xodd + q * xeven;
      s -= p;
      if (s < -1e-10) break;
      if (s >= 0 && Math.abs(2 * s * (xodd - godd)) < 1e-13) break;
    }
  } else tnc = 0;
  tnc += normCdf(-del);
  tnc = Math.min(tnc, 1);
  return negdel ? 1 - tnc : tnc;
}

export const TESTS = ['ind', 'one', 'prop', 'corr'];

/** Cohen's h for two proportions (signed: p1 − p2 direction). */
export const cohenH = (p1, p2) => 2 * Math.asin(Math.sqrt(p1)) - 2 * Math.asin(Math.sqrt(p2));

/**
 * Power of a test. kind: 'ind' (two-sample t, n per group), 'one' (one-sample / paired t, n pairs), 'prop' (two proportions,
 * effect = Cohen's h, n per group), 'corr' (effect = r, n pairs). tails 1 or 2; alpha is the overall (two-tailed) or one-tailed level.
 */
export function power(kind, effect, n, alpha, tails = 2) {
  if (!TESTS.includes(kind)) throw new Error(`kind must be one of ${TESTS.join(', ')} (got ${String(kind)}).`);
  needOpenProbability(alpha, 'alpha');
  if (tails !== 1 && tails !== 2) throw new Error('tails must be 1 or 2.');
  needNumber(effect, 'effect'); needNumber(n, 'n');
  const e = Math.abs(effect);
  const ta = tails === 2 ? alpha / 2 : alpha;
  if (kind === 'ind' || kind === 'one') {
    const df = kind === 'ind' ? 2 * n - 2 : n - 1;
    if (!(df >= 1)) return NaN;
    const ncp = kind === 'ind' ? e * Math.sqrt(n / 2) : e * Math.sqrt(n);
    const tc = tInv(1 - ta, df);
    const upper = 1 - nctCdf(tc, df, ncp);
    return tails === 2 ? upper + nctCdf(-tc, df, ncp) : upper;
  }
  let ncp;
  if (kind === 'prop') ncp = e * Math.sqrt(n / 2);
  else { if (!(n > 3)) return NaN; ncp = Math.atanh(Math.min(e, 0.999999999)) * Math.sqrt(n - 3); }
  const zc = normInv(1 - ta);
  const upper = 1 - normCdf(zc - ncp);
  return tails === 2 ? upper + normCdf(-zc - ncp) : upper;
}

export const minN = (kind) => (kind === 'corr' ? 4 : 2);

/** Smallest integer n (≥ minN) whose power reaches `target`; null when more than 1e8 would be needed. Power is monotone in n. */
export function solveN(kind, effect, target, alpha, tails = 2) {
  if (!TESTS.includes(kind)) throw new Error(`kind must be one of ${TESTS.join(', ')} (got ${String(kind)}).`);
  needOpenProbability(alpha, 'alpha');
  if (tails !== 1 && tails !== 2) throw new Error('tails must be 1 or 2.');
  needNumber(effect, 'effect'); needOpenProbability(target, 'target');
  if (effect === 0) throw new Error('effect cannot be 0: no sample size detects a zero effect.');
  let lo = minN(kind);
  if (power(kind, effect, lo, alpha, tails) >= target) return lo;
  let hi = lo * 2;
  while (!(power(kind, effect, hi, alpha, tails) >= target)) {
    lo = hi; hi *= 2;
    if (hi > 1e8) return null;
  }
  while (hi - lo > 1) {
    const mid = Math.floor((lo + hi) / 2);
    if (power(kind, effect, mid, alpha, tails) >= target) hi = mid; else lo = mid;
  }
  return hi;
}

/** Smallest |effect| (d, h or r) detected with the given power at n; null when unreachable. */
export function solveEffect(kind, n, target, alpha, tails = 2) {
  if (!TESTS.includes(kind)) throw new Error(`kind must be one of ${TESTS.join(', ')} (got ${String(kind)}).`);
  needOpenProbability(alpha, 'alpha');
  if (tails !== 1 && tails !== 2) throw new Error('tails must be 1 or 2.');
  needInteger(n, 'n', 2); needOpenProbability(target, 'target');
  let hi = kind === 'corr' ? 0.999999 : kind === 'prop' ? Math.PI : 20;
  if (!(power(kind, hi, n, alpha, tails) >= target)) return null;
  let lo = 0;
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    if (power(kind, mid, n, alpha, tails) >= target) hi = mid; else lo = mid;
  }
  return hi;
}
