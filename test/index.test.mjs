import test from 'node:test';
import assert from 'node:assert/strict';
import * as api from '../src/index.js';

test('index re-exports the public API', () => {
  for (const name of ['mean', 'sampleSD', 'tSf', 'fSf', 'binomPmf', 'ptukey', 'qtukey', 'power', 'solveN', 'ols', 'oneSampleT', 'twoSampleT']) {
    assert.equal(typeof api[name], 'function', name);
  }
});
