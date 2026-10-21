// Reference values: numpy 1.x / scipy 1.13 (numpy.linalg.lstsq, numpy.linalg.inv, numpy.matmul) computed locally; comments give the call.
import test from 'node:test';
import assert from 'node:assert/strict';
import { matMul, transpose, invUpper, ols } from '../src/matrix.js';

const rel = (a, b, r = 1e-9) => assert.ok(Math.abs(a - b) <= r * Math.max(1, Math.abs(b)), `${a} ≉ ${b}`);

test('matrix: matMul and transpose', () => {
  // numpy.array([[1,2],[3,4],[5,6]]) @ numpy.array([[7,8,9],[10,11,12]]) -> [[27,30,33],[61,68,75],[95,106,117]]
  assert.deepEqual(matMul([[1, 2], [3, 4], [5, 6]], [[7, 8, 9], [10, 11, 12]]), [[27, 30, 33], [61, 68, 75], [95, 106, 117]]);
  assert.deepEqual(transpose([[1, 2, 3], [4, 5, 6]]), [[1, 4], [2, 5], [3, 6]]);
});

test('matrix: invUpper', () => {
  // numpy.linalg.inv([[2,1,3],[0,4,5],[0,0,6]]) -> [[0.5,-0.125,-0.14583333333333331],[0,0.25,-0.20833333333333331],[0,0,0.16666666666666666]]
  const I = invUpper([[2, 1, 3], [0, 4, 5], [0, 0, 6]]);
  rel(I[0][0], 0.5); rel(I[0][1], -0.125); rel(I[0][2], -0.14583333333333331);
  rel(I[1][1], 0.25); rel(I[1][2], -0.20833333333333334); rel(I[2][2], 1 / 6); rel(I[1][0], 0);
});

test('matrix: ols against numpy.linalg.lstsq (two predictors)', () => {
  // X = [1, x1, x2] rows (1,1,2),(1,2,1),(1,3,5),(1,4,3),(1,5,6); y=[3,5,4,9,8]
  // numpy.linalg.lstsq(X, y) -> beta=[1.9666666666666637, 2.3166666666666673, -0.9166666666666653], residual sum of squares 1.1500000000000041;
  // numpy.linalg.inv(X.T @ X) diagonal -> [1.1222222222222218, 0.23888888888888882, 0.1388888888888889]
  const X = [[1, 1, 2], [1, 2, 1], [1, 3, 5], [1, 4, 3], [1, 5, 6]]; const y = [3, 5, 4, 9, 8];
  const f = ols(X, y);
  assert.ok(f.ok);
  rel(f.beta[0], 1.9666666666666637); rel(f.beta[1], 2.3166666666666673); rel(f.beta[2], -0.9166666666666653); rel(f.rss, 1.1500000000000041);
  rel(f.xtxInv[0][0], 1.1222222222222218); rel(f.xtxInv[1][1], 0.23888888888888882); rel(f.xtxInv[2][2], 0.1388888888888889); rel(f.xtxInv[1][2], -0.1388888888888889);
  // normal equations hold: XᵀX β = Xᵀy
  const XtX = matMul(transpose(X), X);
  const Xty = matMul(transpose(X), y.map((v) => [v])).map((r) => r[0]);
  for (let i = 0; i < 3; i++) rel(XtX[i].reduce((a, v, j) => a + v * f.beta[j], 0), Xty[i], 1e-9);
  // (XᵀX)⁻¹ from R: XᵀX · inv = I
  const I = matMul(XtX, f.xtxInv);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) rel(I[i][j], i === j ? 1 : 0, 1e-9);
});

test('matrix: ols exact recovery and residual sum of squares', () => {
  // y = 2 + 3x exactly -> beta [2,3], rss 0
  const f = ols([[1, 0], [1, 1], [1, 2], [1, 3]], [2, 5, 8, 11]);
  rel(f.beta[0], 2); rel(f.beta[1], 3); assert.ok(f.rss < 1e-20);
  // y=[0,1,1] on [1, x=0,1,2]: line a=1/6, b=0.5, rss = 1/6 (numpy.linalg.lstsq residuals -> 0.16666666666666666)
  const g = ols([[1, 0], [1, 1], [1, 2]], [0, 1, 1]);
  rel(g.beta[0], 1 / 6); rel(g.beta[1], 0.5); rel(g.rss, 1 / 6);
});

test('matrix: ols stays accurate with a large offset (years 2000–2011)', () => {
  // statsmodels OLS(val, add_constant(2000..2011)) -> const -548.7405594405811, slope 0.2758741258741368
  const yr = Array.from({ length: 12 }, (_, i) => 2000 + i);
  const val = [3.1, 3.4, 3.3, 3.9, 4.2, 4.1, 4.8, 5.0, 5.3, 5.2, 5.9, 6.1];
  const f = ols(yr.map((v) => [1, v]), val);
  rel(f.beta[0], -548.7405594405811, 1e-8); rel(f.beta[1], 0.2758741258741368, 1e-8);
});

test('matrix: rank deficiency is reported, bad shape throws', () => {
  const f = ols([[1, 1, 2], [1, 2, 4], [1, 3, 6], [1, 4, 8]], [1, 2, 3, 5]);
  assert.equal(f.ok, false);
  assert.equal(f.deficient, 2);
  assert.throws(() => ols([[1, 2, 3]], [1]), /Not enough rows/);
});
