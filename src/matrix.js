// Small dense-matrix helpers for least-squares regression. Pure ESM, no DOM.
// Least squares is solved by Householder QR, not by inverting XᵀX: Golub & Van Loan, "Matrix Computations", 4th ed., §5.2
// (Householder QR) and §5.3 (full-rank least squares); the same approach as R lm() (LINPACK dqrdc2) and numpy.linalg.lstsq.
// The covariance matrix of the estimates is σ²·(XᵀX)⁻¹ = σ²·R⁻¹R⁻ᵀ, with R from the QR factorisation (Montgomery, Peck & Vining,
// "Introduction to Linear Regression Analysis", 5th ed., §3.2).

/** Matrix product of arrays of rows. */
export function matMul(A, B) {
  const m = A.length; const p = B.length; const n = B[0].length;
  const C = Array.from({ length: m }, () => new Array(n).fill(0));
  for (let i = 0; i < m; i++) for (let k = 0; k < p; k++) { const a = A[i][k]; if (a !== 0) for (let j = 0; j < n; j++) C[i][j] += a * B[k][j]; }
  return C;
}

/** Transpose of an array of rows. */
export function transpose(A) {
  return A[0].map((_, j) => A.map((row) => row[j]));
}

/** Inverse of an upper-triangular matrix (array of rows) by back substitution. */
export function invUpper(R) {
  const p = R.length;
  const inv = Array.from({ length: p }, () => new Array(p).fill(0));
  for (let j = 0; j < p; j++) {
    inv[j][j] = 1 / R[j][j];
    for (let i = j - 1; i >= 0; i--) {
      let s = 0;
      for (let k = i + 1; k <= j; k++) s += R[i][k] * inv[k][j];
      inv[i][j] = -s / R[i][i];
    }
  }
  return inv;
}

/** Relative size below which a diagonal element of R is treated as zero (column is a linear combination of earlier ones). */
export const RANK_TOL = 1e-10;

/**
 * Ordinary least squares y ≈ X·β by Householder QR.
 * X: array of m rows with p columns (include the column of ones yourself), y: array of m numbers, m ≥ p.
 * Returns { ok: true, beta, rss, xtxInv, rank: p } or { ok: false, deficient: j } where j (0-based) is the first column that is
 * (numerically) a linear combination of the columns before it.
 */
export function ols(X, y) {
  if (!Array.isArray(X) || X.length === 0 || !Array.isArray(X[0]) || X[0].length === 0) throw new Error('X must be a non-empty array of rows.');
  if (!Array.isArray(y) || y.length !== X.length) throw new Error('y must be an array with one value per row of X.');
  const m = X.length; const p = X[0].length;
  for (let i = 0; i < m; i++) {
    if (!Array.isArray(X[i]) || X[i].length !== p) throw new Error(`Row ${i} of X must have ${p} columns.`);
    for (let j = 0; j < p; j++) if (!Number.isFinite(X[i][j])) throw new Error(`X[${i}][${j}] is not a finite number.`);
    if (!Number.isFinite(y[i])) throw new Error(`y[${i}] is not a finite number.`);
  }
  if (m < p) throw new Error(`Not enough rows (${m}) for ${p} unknowns.`);
  const cols = [];
  const norm0 = [];
  for (let j = 0; j < p; j++) {
    const c = new Float64Array(m);
    let s = 0;
    for (let i = 0; i < m; i++) { c[i] = X[i][j]; s += c[i] * c[i]; }
    cols.push(c); norm0.push(Math.sqrt(s));
  }
  const b = Float64Array.from(y);
  for (let k = 0; k < p; k++) {
    let nrm = 0;
    for (let i = k; i < m; i++) nrm += cols[k][i] * cols[k][i];
    nrm = Math.sqrt(nrm);
    if (!(nrm > RANK_TOL * norm0[k]) || nrm === 0) return { ok: false, deficient: k };
    const alpha = cols[k][k] > 0 ? -nrm : nrm;
    const v = new Float64Array(m);
    for (let i = k; i < m; i++) v[i] = cols[k][i];
    v[k] -= alpha;
    let vv = 0;
    for (let i = k; i < m; i++) vv += v[i] * v[i];
    if (vv > 0) {
      const apply = (c) => {
        let d = 0;
        for (let i = k; i < m; i++) d += v[i] * c[i];
        const f = (2 * d) / vv;
        for (let i = k; i < m; i++) c[i] -= f * v[i];
      };
      for (let j = k + 1; j < p; j++) apply(cols[j]);
      apply(b);
    }
    cols[k][k] = alpha;
    for (let i = k + 1; i < m; i++) cols[k][i] = 0;
  }
  const R = Array.from({ length: p }, (_, i) => Array.from({ length: p }, (_, j) => (j >= i ? cols[j][i] : 0)));
  const beta = new Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) {
    let s = b[i];
    for (let j = i + 1; j < p; j++) s -= R[i][j] * beta[j];
    beta[i] = s / R[i][i];
  }
  let rss = 0;
  for (let i = p; i < m; i++) rss += b[i] * b[i];
  const Ri = invUpper(R);
  const xtxInv = matMul(Ri, transpose(Ri));
  return { ok: true, beta, rss, xtxInv, rank: p };
}
