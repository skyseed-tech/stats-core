// Reference values: scipy 1.13.1 scipy.stats.studentized_range.cdf(q, k, df) / .ppf(p, k, df) (same quantity as R ptukey / qtukey).
// ptukey is a port of R's Gauss–Legendre algorithm (accuracy about 1e-8), so tolerances are 1e-7 absolute (1e-6 for quantiles).
import test from 'node:test';
import assert from 'node:assert/strict';
import { ptukey, ptukeySf, qtukey } from '../src/tukey.js';

const close = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol, `${a} ≉ ${b}`);

test('ptukey matches scipy studentized_range.cdf', () => {
  // [q, k, df, scipy value]
  const cases = [
    [3.5, 3, 20, 0.94410815091825], [3.877, 3, 20, 0.9677197752303514], [1, 2, 5, 0.48891591956971947],
    [5, 4, 10, 0.9765540036340024], [2.0, 10, 30, 0.0868471908954611], [4.5, 5, 100, 0.9836657545708238],
    [3.314, 3, 10, 0.9047960349430152], [7, 8, 3, 0.9069165647373614], [3.0, 3, 100000, 0.9144574283450421],
    [0.5, 3, 12, 0.0662089800462016],
  ];
  for (const [q, k, df, want] of cases) close(ptukey(q, k, df), want, 1e-7);
});

test('qtukey matches scipy studentized_range.ppf', () => {
  const cases = [
    [0.95, 3, 20, 3.577934725220135], [0.95, 2, 10, 3.1510641833294026], [0.99, 5, 12, 5.836308366059154],
    [0.95, 10, 30, 4.824141286183029], [0.9, 4, 60, 3.311904339491206], [0.95, 3, 2000, 3.31696747982079],
  ];
  for (const [p, k, df, want] of cases) close(qtukey(p, k, df), want, 1e-6);
});

test('qtukey inverts ptukey; k = 2 equals sqrt(2) × the t quantile (Tukey = t test)', () => {
  for (const [p, k, df] of [[0.5, 4, 15], [0.975, 6, 40], [0.05, 3, 8]]) close(ptukey(qtukey(p, k, df), k, df), p, 1e-9);
  // scipy: studentized_range.ppf(0.95, 2, 10) = 3.1510641833 = sqrt(2) * t.ppf(0.975, 10) = 1.41421356 * 2.22813886
  close(qtukey(0.95, 2, 10), Math.SQRT2 * 2.2281388519649385, 1e-6);
});

test('edge cases', () => {
  assert.equal(ptukey(0, 3, 10), 0);
  assert.equal(ptukey(-1, 3, 10), 0);
  assert.equal(ptukey(Infinity, 3, 10), 1);
  assert.ok(Number.isNaN(ptukey(3, 3, 1)));      // df must be at least 2
  assert.ok(Number.isNaN(ptukey(3, 1, 10)));     // needs at least 2 means
  assert.equal(qtukey(0, 3, 10), 0);
  assert.equal(qtukey(1, 3, 10), Infinity);
  assert.ok(Number.isNaN(qtukey(1.5, 3, 10)));
  close(ptukeySf(3.5, 3, 20), 1 - 0.94410815091825, 1e-7);
  assert.ok(ptukeySf(40, 3, 20) >= 0);
});
