import test from 'node:test';
import assert from 'node:assert/strict';
import { summarize, oneSampleT, twoSampleT, tPValue, tCritical, tCritTwo, binomPmf, binomTestTwoSided, poisPmf, power, solveN, solveEffect, ols } from '../src/index.js';

test('summarize: rejects non-arrays, short samples and non-finite values', () => {
  assert.throws(() => summarize('1 2 3'), /must be an array/);
  assert.throws(() => summarize([]), /at least 2 values \(got 0\)/);
  assert.throws(() => summarize([1]), /at least 2 values/);
  assert.throws(() => summarize([1, NaN, 3]), /xs\[1\] is not a finite number/);
  assert.throws(() => summarize([1, '2']), /xs\[1\]/);
  assert.deepEqual(summarize([1, 2, 3]), { n: 3, mean: 2, sd: 1 });
});

test('oneSampleT: validates inputs', () => {
  const ok = { mean: 5, sd: 2, n: 10 };
  assert.throws(() => oneSampleT({ ...ok, n: 1 }), /n must be a whole number of at least 2/);
  assert.throws(() => oneSampleT({ ...ok, n: 2.5 }), /whole number/);
  assert.throws(() => oneSampleT({ ...ok, sd: -1 }), /sd cannot be negative/);
  assert.throws(() => oneSampleT({ ...ok, sd: 0 }), /identical/);
  assert.throws(() => oneSampleT({ ...ok, mean: NaN }), /mean must be a finite number/);
  assert.throws(() => oneSampleT({ ...ok, tail: 'both' }), /tail must be/);
  assert.throws(() => oneSampleT({ ...ok, alpha: 1.5 }), /alpha must be strictly between 0 and 1/);
  assert.ok(Number.isFinite(oneSampleT(ok).p));
});

test('twoSampleT: validates inputs', () => {
  const ok = { m1: 5, s1: 1, n1: 10, m2: 4, s2: 1, n2: 10 };
  assert.throws(() => twoSampleT({ ...ok, n2: 1 }), /n2 must be a whole number/);
  assert.throws(() => twoSampleT({ ...ok, s1: -1 }), /s1 cannot be negative/);
  assert.throws(() => twoSampleT({ ...ok, s1: 0, s2: 0 }), /both 0/);
  assert.throws(() => twoSampleT({ ...ok, m2: Infinity }), /m2 must be a finite number/);
  assert.throws(() => twoSampleT({ ...ok, tail: 'x' }), /tail must be/);
  assert.ok(Number.isFinite(twoSampleT({ ...ok, s1: 0 }).p));
});

test('p-values and critical values: validate inputs', () => {
  assert.throws(() => tPValue(NaN, 5), /stat must be a finite number/);
  assert.throws(() => tPValue(1, 0), /df must be greater than 0/);
  assert.throws(() => tPValue(1, 5, 'up'), /tail must be/);
  assert.throws(() => tCritTwo(0, 10), /alpha/);
  assert.throws(() => tCritical(0.05, -1, 'two'), /df/);
});

test('discrete distributions: validate parameters', () => {
  assert.throws(() => binomPmf(1, 10.5, 0.5), /n must be a whole number/);
  assert.throws(() => binomPmf(1, 10, 1.5), /p must be between 0 and 1/);
  assert.throws(() => binomPmf(NaN, 10, 0.5), /x must be a finite number/);
  assert.throws(() => binomTestTwoSided(1, -1, 0.5), /n must be/);
  assert.throws(() => poisPmf(1, -2), /lambda cannot be negative/);
  assert.throws(() => poisPmf(1, NaN), /lambda/);
});

test('power and sample size: validate inputs', () => {
  assert.throws(() => power('anova', 0.5, 20, 0.05), /kind must be one of/);
  assert.throws(() => power('ind', 0.5, 20, 0), /alpha/);
  assert.throws(() => power('ind', 0.5, 20, 0.05, 3), /tails must be 1 or 2/);
  assert.throws(() => power('ind', NaN, 20, 0.05), /effect/);
  assert.throws(() => solveN('ind', 0, 0.8, 0.05), /effect cannot be 0/);
  assert.throws(() => solveN('ind', 0.5, 1, 0.05), /target/);
  assert.throws(() => solveEffect('one', 1, 0.8, 0.05), /n must be a whole number of at least 2/);
});

test('ols: validates shape and values', () => {
  assert.throws(() => ols([], []), /non-empty array of rows/);
  assert.throws(() => ols('x', [1]), /non-empty array of rows/);
  assert.throws(() => ols([[1, 2], [1, 3]], [1]), /one value per row/);
  assert.throws(() => ols([[1, 2], [1]], [1, 2]), /Row 1 of X must have 2 columns/);
  assert.throws(() => ols([[1, NaN], [1, 3]], [1, 2]), /X\[0\]\[1\]/);
  assert.throws(() => ols([[1, 2], [1, 3]], [1, NaN]), /y\[1\]/);
  assert.throws(() => ols([[1, 2, 3]], [1]), /Not enough rows/);
});
