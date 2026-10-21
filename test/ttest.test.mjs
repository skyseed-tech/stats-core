// Sources: scipy 1.13 (run locally) — stats.ttest_1samp, ttest_ind(equal_var=True/False), ttest_ind_from_stats, stats.t.ppf/sf.
import test from 'node:test';
import assert from 'node:assert/strict';
import { oneSampleT, twoSampleT, tPValue, tCritical, tCritTwo, effectKey, summarize } from '../src/ttest.js';

const rel = (a, b, r = 1e-6) => assert.ok(Math.abs(a - b) <= r * Math.abs(b), `${a} ≉ ${b}`);

test('ttc: one-sample from raw data matches scipy', () => {
  // scipy 1.13: stats.ttest_1samp([102,98,105,110,99,101,97,104,108,100], 100) -> t=1.76505, p=0.111373, df=9, CI=(99.32406, 105.47594)
  const s = summarize([102, 98, 105, 110, 99, 101, 97, 104, 108, 100]);
  const r = oneSampleT({ ...s, mu0: 100, tail: 'two', alpha: 0.05 });
  rel(r.t, 1.7650452162436605); rel(r.p, 0.11137324638262684); assert.equal(r.df, 9);
  rel(r.ciLow, 99.32405772901154); rel(r.ciHigh, 105.47594227098847);
  rel(r.d, 2.4 / 4.2998707990925595);
});

test('ttc: one-sample one-sided', () => {
  // scipy 1.13: stats.ttest_1samp(x, 98, alternative='greater') -> t=3.23592, p=0.00511264962
  const s = summarize([102, 98, 105, 110, 99, 101, 97, 104, 108, 100]);
  const r = oneSampleT({ ...s, mu0: 98, tail: 'greater' });
  rel(r.t, 3.235916229780041); rel(r.p, 0.00511264962005501, 1e-5);
});

test('ttc: Welch and pooled from raw data (Wikipedia / scipy example data)', () => {
  const a1 = [27.5, 21.0, 19.0, 23.6, 17.0, 17.9, 16.9, 20.1, 21.9, 22.6, 23.1, 19.6, 19.0, 21.7, 21.4];
  const a2 = [27.1, 22.0, 20.8, 23.4, 23.4, 23.5, 25.8, 22.0, 24.8, 20.2, 21.9, 22.1, 22.9, 20.5, 24.4];
  const x = summarize(a1); const y = summarize(a2);
  const base = { m1: x.mean, s1: x.sd, n1: x.n, m2: y.mean, s2: y.sd, n2: y.n };
  // scipy 1.13: stats.ttest_ind(a1,a2,equal_var=False) -> t=-2.455356, p=0.021378, df=24.98853, CI=(-3.984096, -0.349237)
  const w = twoSampleT({ ...base, pooled: false });
  rel(w.t, -2.455356398286006); rel(w.p, 0.021378001462866985); rel(w.df, 24.988529290231416);
  rel(w.ciLow, -3.9840962671409264); rel(w.ciHigh, -0.34923706619240935);
  // scipy 1.13: stats.ttest_ind(a1,a2,equal_var=True) -> t=-2.455356, p=0.0205445, df=28, CI=(-3.974231, -0.359102)
  const p = twoSampleT({ ...base, pooled: true });
  rel(p.p, 0.020544522734125933); assert.equal(p.df, 28); rel(p.ciLow, -3.9742313350740126); rel(p.ciHigh, -0.3591019982593229);
  // alternative='less' (Welch) -> p=0.010689 ; alternative='greater' (pooled) -> p=0.98972774
  rel(twoSampleT({ ...base, tail: 'less' }).p, 0.010689000731433492);
  rel(twoSampleT({ ...base, pooled: true, tail: 'greater' }).p, 0.9897277386329371);
});

test('ttc: two-sample from summary statistics', () => {
  // scipy 1.13: stats.ttest_ind_from_stats(5.2,1.1,20,4.6,1.3,22,equal_var=True) -> t=1.606125, p=0.116113
  const r = twoSampleT({ m1: 5.2, s1: 1.1, n1: 20, m2: 4.6, s2: 1.3, n2: 22, pooled: true });
  rel(r.t, 1.606124685622228); rel(r.p, 0.1161133432870854); rel(r.sp, 1.2091319200153472);
  rel(r.d, 0.496223770184137);
  // equal_var=False -> t=1.619151, p=0.113311 ; df by the Welch–Satterthwaite formula = 39.81089
  const w = twoSampleT({ m1: 5.2, s1: 1.1, n1: 20, m2: 4.6, s2: 1.3, n2: 22 });
  rel(w.t, 1.6191505367860586); rel(w.p, 0.11331103069921064); rel(w.df, 39.81089394959134);
  // delta0 shifts the statistic: (0.6-0.5)/se
  rel(twoSampleT({ m1: 5.2, s1: 1.1, n1: 20, m2: 4.6, s2: 1.3, n2: 22, delta0: 0.5, pooled: true }).t, 0.1 / 0.3735700007422243);
});

test('ttc: critical values and p-values', () => {
  // scipy 1.13: stats.t.ppf(0.975, 10) = 2.2281388519649385; t.ppf(0.975, 2.5) = 3.5746548420118778; t.ppf(.95, 15.5) = 1.7493443916588804
  rel(tCritTwo(0.05, 10), 2.2281388519649385, 1e-7);
  rel(tCritTwo(0.05, 2.5), 3.5746548420118778, 1e-6);
  rel(tCritical(0.05, 15.5, 'greater').high, 1.7493443916588804, 1e-6);
  assert.equal(tCritical(0.05, 10, 'less').high, null);
  // scipy: stats.t.sf(2.2694, 15.5) = 0.018954964292269776
  rel(tPValue(2.2694, 15.5, 'greater'), 0.018954964292269776, 1e-7);
  rel(tPValue(-2.2694, 15.5, 'two'), 2 * 0.018954964292269776, 1e-7);
});

test('ttc: effect-size labels', () => {
  assert.equal(effectKey(0.1), 'none'); assert.equal(effectKey(-0.3), 'small'); assert.equal(effectKey(0.6), 'medium'); assert.equal(effectKey(0.9), 'large');
});
