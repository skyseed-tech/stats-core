// Sources: normal CDF/inverse: A&S 7.1 and Acklam (checked against 1.95996 / 2.32635 etc.); d2: NIST/Montgomery; F distribution checked against scipy.
import test from 'node:test';
import assert from 'node:assert/strict';
import { mean, sampleSD, popSD, normCdf, normSf, normInv, normInvSf, fSf, lnGamma, parseNumberList, parseOne, D2 } from '../src/stats.js';

const close = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) <= tol, `${a} ≉ ${b}`);
const rel = (a, b, r = 1e-8) => assert.ok(Math.abs(a - b) <= r * Math.abs(b), `${a} ≉ ${b}`);

test('mean and standard deviations (NIST-style hand values)', () => {
  // 2,4,4,4,5,5,7,9: mean 5, population SD 2 (classic example), sample SD = sqrt(32/7)
  const xs = [2, 4, 4, 4, 5, 5, 7, 9];
  close(mean(xs), 5);
  close(popSD(xs), 2);
  close(sampleSD(xs), Math.sqrt(32 / 7));
  assert.ok(Number.isNaN(sampleSD([3])));
});

test('Φ against published table values (A&S 26.1 / standard normal tables)', () => {
  close(normCdf(0), 0.5, 1e-15);
  close(normCdf(1), 0.8413447460685429, 1e-12);
  close(normCdf(1.96), 0.9750021048517795, 1e-12);
  close(normCdf(-1.5), 0.0668072012688581, 1e-12);
  close(normCdf(3), 0.9986501019683699, 1e-12);
  rel(normSf(4.5), 3.3976731247300535e-6, 1e-9);
  rel(normSf(6), 9.865876450376946e-10, 1e-9);
  rel(normCdf(-8), 6.220960574271785e-16, 1e-8);
  assert.equal(normCdf(Infinity), 1);
});

test('Φ⁻¹ round trip and known quantiles', () => {
  close(normInv(0.975), 1.959963984540054, 1e-12);
  close(normInv(0.5), 0, 1e-15);
  close(normInv(0.8413447460685429), 1, 1e-10);
  rel(normInvSf(3.4e-6), 4.4999, 1e-4);
  for (const p of [1e-9, 1e-6, 0.001, 0.02, 0.3, 0.5, 0.9, 0.99, 0.999999]) rel(normCdf(normInv(p)), p, 1e-9);
  assert.equal(normInv(1), Infinity);
  assert.ok(Number.isNaN(normInv(1.2)));
});

test('F distribution upper tail (checked against scipy.stats.f.sf)', () => {
  close(fSf(1, 1, 1), 0.5, 1e-12);
  close(fSf(3, 2, 2), 1 / 4, 1e-12); // d1 = d2 = 2: 1/(1+F)
  rel(fSf(4.0, 5, 20), 0.011183751855265593, 1e-9); // scipy.stats.f.sf(4, 5, 20)
  close(lnGamma(5), Math.log(24), 1e-12);
  close(lnGamma(0.5), Math.log(Math.sqrt(Math.PI)), 1e-12);
});

test('parseNumberList and parseOne', () => {
  assert.deepEqual(parseNumberList('1, 2.5  3\n4;5\t−6 1e2 .5').values, [1, 2.5, 3, 4, 5, -6, 100, 0.5]);
  assert.deepEqual(parseNumberList('1 x 2 3,,').bad, ['x']);
  assert.ok(Number.isNaN(parseOne('')));
  assert.equal(parseOne(' 1,000 '), 1000);
});

test('d2 constants for subgroup sizes 2-10', () => {
  assert.deepEqual(Object.keys(D2).map(Number), [2, 3, 4, 5, 6, 7, 8, 9, 10]);
  assert.equal(D2[5], 2.326);
});
