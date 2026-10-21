// Student / Welch t-test maths from summary statistics. Pure ESM, no DOM.
// Sources:
//  * One-sample t test: NIST/SEMATECH e-Handbook of Statistical Methods 1.3.5.2 and 7.2.2 — t = (x̄ − μ0)/(s/√n), df = n − 1
//    (https://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm). Same as R t.test(x, mu=) / scipy.stats.ttest_1samp.
//  * Two-sample test: NIST 7.3.1 "Do two samples have the same mean?" — Welch–Satterthwaite df = (s1²/n1 + s2²/n2)² /
//    [(s1²/n1)²/(n1−1) + (s2²/n2)²/(n2−1)] (non-integer; same as R t.test(var.equal=FALSE) / scipy.stats.ttest_ind(equal_var=False)),
//    pooled (Student) variant sp² = [(n1−1)s1² + (n2−1)s2²]/(n1+n2−2), df = n1 + n2 − 2 (R var.equal=TRUE / equal_var=True).
//  * Confidence interval for the (difference of) mean(s): estimate ± t(1−α/2, df)·SE (two-sided, always, whatever the tail of the test).
//  * Effect size: Cohen (1988), Statistical Power Analysis for the Behavioral Sciences, 2nd ed., ch. 2: d = (x̄1 − x̄2)/s_pooled,
//    one-sample d = (x̄ − μ0)/s, paired d_z = d̄/s_d; benchmarks 0.2 small, 0.5 medium, 0.8 large.
import { mean, sampleSD, tCdf, tSf, tInv } from './stats.js';
import { needNumber, needPositive, needNonNegative, needInteger, needOpenProbability, needTail, needSample } from './validate.js';

export const TAILS = ['two', 'greater', 'less'];

/** p-value of a t statistic: tail 'two' = 2·P(T ≥ |t|), 'greater' = P(T ≥ t), 'less' = P(T ≤ t). */
export function tPValue(stat, df, tail = 'two') {
  needNumber(stat, 'stat'); needPositive(df, 'df'); needTail(tail);
  if (tail === 'greater') return tSf(stat, df);
  if (tail === 'less') return tCdf(stat, df);
  return Math.min(1, 2 * tSf(Math.abs(stat), df));
}
/** p-value of a z statistic with the same tail names. */
export function zPValue(stat, normSfFn, tail = 'two') {
  needNumber(stat, 'stat'); needTail(tail);
  if (typeof normSfFn !== 'function') throw new Error('normSfFn must be a function.');
  if (tail === 'greater') return normSfFn(stat);
  if (tail === 'less') return normSfFn(-stat);
  return Math.min(1, 2 * normSfFn(Math.abs(stat)));
}
/** Two-sided t critical value t(1−α/2, df). */
export function tCritTwo(alpha, df) {
  needOpenProbability(alpha, 'alpha'); needPositive(df, 'df');
  return -tInv(alpha / 2, df);
}
/** Critical value(s) of a t test at α for the tail: { low, high } (null when not used). */
export function tCritical(alpha, df, tail = 'two') {
  needOpenProbability(alpha, 'alpha'); needPositive(df, 'df'); needTail(tail);
  if (tail === 'greater') return { low: null, high: -tInv(alpha, df) };
  if (tail === 'less') return { low: tInv(alpha, df), high: null };
  const c = tCritTwo(alpha, df);
  return { low: -c, high: c };
}

/** Cohen's benchmarks → key of the size label. */
export function effectKey(d) {
  const a = Math.abs(d);
  return a < 0.2 ? 'none' : a < 0.5 ? 'small' : a < 0.8 ? 'medium' : 'large';
}

/** Mean, sample SD and n of a list. */
export function summarize(xs) {
  needSample(xs, 'xs', 2);
  return { n: xs.length, mean: mean(xs), sd: sampleSD(xs) };
}

/**
 * One-sample t test from summary statistics (also the paired test, applied to the differences).
 * { mean, sd, n, mu0, tail, alpha } → { n, mean, sd, se, df, t, p, diff, ciLow, ciHigh, d }
 */
export function oneSampleT({ mean: m, sd, n, mu0 = 0, tail = 'two', alpha = 0.05 }) {
  needNumber(m, 'mean'); needNonNegative(sd, 'sd'); needInteger(n, 'n', 2); needNumber(mu0, 'mu0'); needTail(tail); needOpenProbability(alpha, 'alpha');
  if (sd === 0) throw new Error('sd is 0: all values are identical, so the t statistic is undefined.');
  const se = sd / Math.sqrt(n);
  const df = n - 1;
  const stat = (m - mu0) / se;
  const c = tCritTwo(alpha, df);
  return {
    n, mean: m, sd, se, df, t: stat, p: tPValue(stat, df, tail), diff: m - mu0,
    ciLow: m - c * se, ciHigh: m + c * se, tCrit: c, d: (m - mu0) / sd,
  };
}

/**
 * Two-sample t test from summary statistics.
 * { m1, s1, n1, m2, s2, n2, delta0, pooled, tail, alpha } → { se, df, t, p, diff, ciLow, ciHigh, sp, d, … }
 * The tested quantity is μ1 − μ2 against delta0. Cohen's d always uses the pooled SD sp.
 */
export function twoSampleT({ m1, s1, n1, m2, s2, n2, delta0 = 0, pooled = false, tail = 'two', alpha = 0.05 }) {
  needNumber(m1, 'm1'); needNumber(m2, 'm2'); needNonNegative(s1, 's1'); needNonNegative(s2, 's2');
  needInteger(n1, 'n1', 2); needInteger(n2, 'n2', 2); needNumber(delta0, 'delta0'); needTail(tail); needOpenProbability(alpha, 'alpha');
  if (s1 === 0 && s2 === 0) throw new Error('s1 and s2 are both 0: the t statistic is undefined.');
  const v1 = s1 * s1 / n1; const v2 = s2 * s2 / n2;
  const sp2 = ((n1 - 1) * s1 * s1 + (n2 - 1) * s2 * s2) / (n1 + n2 - 2);
  const sp = Math.sqrt(sp2);
  let se; let df;
  if (pooled) { se = sp * Math.sqrt(1 / n1 + 1 / n2); df = n1 + n2 - 2; } else {
    se = Math.sqrt(v1 + v2);
    df = (v1 + v2) ** 2 / (v1 * v1 / (n1 - 1) + v2 * v2 / (n2 - 1));
  }
  const diff = m1 - m2;
  const stat = (diff - delta0) / se;
  const c = tCritTwo(alpha, df);
  return {
    se, df, t: stat, p: tPValue(stat, df, tail), diff, ciLow: diff - c * se, ciHigh: diff + c * se, tCrit: c,
    sp, d: diff / sp, v1, v2,
  };
}
