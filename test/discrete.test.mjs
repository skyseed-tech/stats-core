// Reference values: scipy 1.13 (run locally), scipy.stats.binom / poisson / binomtest. Tolerance 1e-10 relative (the algorithm is Loader's, ~1e-14).
import test from 'node:test';
import assert from 'node:assert/strict';
import { binomPmf, binomCdf, binomSf, binomRange, binomTestTwoSided, poisPmf, poisCdf, poisSf, poisRange, eventBounds, normalApprox, stirlerr } from '../src/discrete.js';

const rel = (a, b, r = 1e-10) => assert.ok(Math.abs(a - b) <= r * Math.abs(b), `${a} ≉ ${b}`);

test('discrete: binomial pmf / cdf / sf', () => {
  rel(binomPmf(3, 10, 0.5), 0.11718750000000004);          // stats.binom.pmf(3,10,.5)
  rel(binomCdf(3, 10, 0.5), 0.171875);                      // stats.binom.cdf(3,10,.5)
  rel(binomSf(2, 10, 0.5), 0.9453125);                      // stats.binom.sf(2,10,.5)
  rel(binomRange(4, 8, 20, 0.3), 0.7795817326192909);       // cdf(8,20,.3)-cdf(3,20,.3)
  rel(binomSf(59, 100, 0.5), 0.028443966820490392);         // stats.binom.sf(59,100,.5)
  rel(binomPmf(0, 100, 0.1), 2.656139888758749e-05);        // stats.binom.pmf(0,100,.1)
  rel(binomCdf(5, 100, 0.5), 6.26162256269268e-23);         // stats.binom.cdf(5,100,.5)
});

test('discrete: binomial with n = 10^6', () => {
  rel(binomPmf(500000, 1e6, 0.5), 0.0007978843613317499, 1e-12);   // stats.binom.pmf(500000,10**6,.5)
  rel(binomSf(1099, 1e6, 0.001), 0.0009574669700286353, 1e-9);     // stats.binom.sf(1099,10**6,.001)
  rel(binomCdf(990, 1e6, 0.001), 0.38370193487995213, 1e-9);       // stats.binom.cdf(990,10**6,.001)
});

test('discrete: exact binomial test', () => {
  rel(binomTestTwoSided(7, 10, 0.5), 0.34375);                     // stats.binomtest(7,10,.5).pvalue
  rel(binomTestTwoSided(60, 100, 0.5), 0.056887933640980784);      // stats.binomtest(60,100,.5).pvalue
  rel(binomTestTwoSided(3, 20, 0.3), 0.22041826738070913);         // stats.binomtest(3,20,.3).pvalue
  rel(binomTestTwoSided(530, 1000, 0.5), 0.06202319509836341);     // stats.binomtest(530,1000,.5).pvalue
  assert.ok(binomTestTwoSided(2, 10, 0.2) <= 1 && binomTestTwoSided(2, 10, 0.2) > 0.999999); // binomtest(2,10,.2).pvalue = 1.0
});

test('discrete: Poisson pmf / cdf / sf', () => {
  rel(poisPmf(2, 3.5), 0.18495897346170082);               // stats.poisson.pmf(2,3.5)
  rel(poisCdf(2, 3.5), 0.3208471988621341);                 // stats.poisson.cdf(2,3.5)
  rel(poisSf(7, 12), 0.9104955031598242);                   // stats.poisson.sf(7,12)
  rel(poisRange(9, 15, 12), 0.6893878706827202);            // cdf(15,12)-cdf(8,12)
  rel(poisSf(1099, 1000), 0.0009626304058665588, 1e-9);    // stats.poisson.sf(1099,1000)
  rel(poisPmf(1000, 1000), 0.01261461134870819, 1e-11);    // stats.poisson.pmf(1000,1000)
  rel(poisSf(40, 10), 1.777341749349938e-13, 1e-8);        // stats.poisson.sf(40,10)
  rel(poisPmf(100, 100.5), 0.03981136695234592);           // stats.poisson.pmf(100,100.5)
});

test('discrete: helpers', () => {
  assert.deepEqual(eventBounds('lt', 5, 5, 10), [0, 4]);
  assert.deepEqual(eventBounds('gt', 5, 5, 10), [6, 10]);
  assert.deepEqual(eventBounds('between', 2, 7, 10), [2, 7]);
  assert.equal(normalApprox(3, 2, 5, 1), null);
  assert.ok(Math.abs(normalApprox(0, Infinity, 5, 1) - 1) < 1e-4);
  assert.ok(Math.abs(stirlerr(20) - 1 / 240.5) < 1e-4); // ≈ 1/(12 n)
  assert.equal(binomPmf(5, 3, 0.5), 0);
  assert.equal(poisPmf(0, 0), 1);
  assert.equal(binomRange(0, 5, 10, 0), 1);
});
