// Noncentral t CDF (AS 243) against scipy 1.13 stats.nct.cdf(x, df, nc); power helpers against statsmodels 0.14.
import test from 'node:test';
import assert from 'node:assert/strict';
import { nctCdf, power, solveN, solveEffect, cohenH } from '../src/power.js';

const rel = (a, b, r = 1e-8) => assert.ok(Math.abs(a - b) <= r * Math.abs(b), `${a} ≉ ${b}`);

test('power.js: noncentral t CDF equals scipy.stats.nct.cdf', () => {
  // scipy 1.13: nct.cdf(x, df, nc)
  rel(nctCdf(2.0, 18, 2.5), 0.3080015288017259, 1e-9);
  rel(nctCdf(-2.0, 18, 2.5), 1.0166811473477644e-05, 1e-6);
  rel(nctCdf(2.1009, 18, 3.1623), 0.15154285052988187, 1e-9);
  rel(nctCdf(-2.1009, 18, 3.1623), 3.3827885248527423e-07, 1e-5);
  rel(nctCdf(1.7, 5, 1), 0.7064946909056553, 1e-9);
  rel(nctCdf(0.3, 40, -1), 0.9027543304347568, 1e-9);
  rel(nctCdf(5, 100, 8), 0.0022733227655674522, 1e-9);
  rel(nctCdf(2.0, 2, 0.5), 0.8242224279725603, 1e-9);
  rel(nctCdf(-2.0, 98, 2.8), 1.0294160291968169e-06, 1e-5);
  rel(nctCdf(1.96, 1000, 3), 0.14928791553444706, 1e-9);
});

test('power.js: ncp = 0 is the central t; edge arguments', () => {
  rel(nctCdf(1.5, 10, 0), 0.9177463367772799, 1e-10); // scipy t.cdf(1.5, 10)
  assert.equal(nctCdf(Infinity, 5, 1), 1);
  assert.equal(nctCdf(-Infinity, 5, 1), 0);
  assert.ok(Number.isNaN(nctCdf(1, -1, 1)));
});

test('power.js: t-test power, sample size and minimal effect equal statsmodels', () => {
  // TTestIndPower().power(0.5, 64, 0.05) = 0.8014595500498423 ; (0.5, 30, .05) = 0.4778965195865288
  rel(power('ind', 0.5, 64, 0.05, 2), 0.8014595500498423, 1e-8);
  rel(power('ind', 0.5, 30, 0.05, 2), 0.4778965195865288, 1e-8);
  // alternative='larger': power(0.5, 20, 0.05) = 0.46337434869922856
  rel(power('ind', 0.5, 20, 0.05, 1), 0.46337434869922856, 1e-8);
  // TTestPower().power(0.3, 50, 0.01) = 0.299469489332798
  rel(power('one', 0.3, 50, 0.01, 2), 0.299469489332798, 1e-8);
  // solve_power(0.5, None, 0.05, 0.8) = 63.7656 -> 64 ; (0.2) -> 393.41 -> 394 ; TTestPower(0.8, 0.95) = 22.32 -> 23
  assert.equal(solveN('ind', 0.5, 0.8, 0.05, 2), 64);
  assert.equal(solveN('ind', 0.2, 0.8, 0.05, 2), 394);
  assert.equal(solveN('one', 0.8, 0.95, 0.05, 2), 23);
  // TTestIndPower d=0.3 one-sided power 0.8 -> 138.07 -> 139
  assert.equal(solveN('ind', 0.3, 0.8, 0.05, 1), 139);
  // minimal effect: TTestIndPower().solve_power(None, 30, 0.05, 0.8) = 0.7356198 (statsmodels' root finder is good to ~1e-6)
  rel(solveEffect('ind', 30, 0.8, 0.05, 2), 0.7356198427826622, 1e-5);
});

test('power.js: proportions (Cohen h) and correlation', () => {
  // proportion_effectsize(0.6, 0.5) = 0.20135792079033066
  rel(Math.abs(cohenH(0.6, 0.5)), 0.20135792079033066, 1e-12);
  // NormalIndPower().power(h, 200, 0.05) = 0.5214145419211704
  rel(power('prop', 0.20135792079033066, 200, 0.05, 2), 0.5214145419211704, 1e-9);
  assert.equal(solveN('prop', 0.20135792079033066, 0.8, 0.05, 2), 388);
  // Fisher z: n = ((z_.975 + z_.8)/atanh(r))² + 3 = 84.9278 -> 85 ; power(r=.3, n=50) = 0.5643676389635401 (scipy norm)
  assert.equal(solveN('corr', 0.3, 0.8, 0.05, 2), 85);
  rel(power('corr', 0.3, 50, 0.05, 2), 0.5643676389635401, 1e-9);
});
