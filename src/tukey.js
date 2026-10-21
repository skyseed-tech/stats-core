// Studentized range distribution (Tukey HSD): ptukey / qtukey. Pure maths, no DOM.
// Source: a port of R's ptukey.c (R Core Team, "stats" package), which implements the algorithm of
//   Copenhaver & Holland (1988), "Multiple comparisons of simple effects in the two-way analysis of variance model",
//   Computational Statistics & Data Analysis 6, 295–302 (AS 190 lineage, Gauss–Legendre quadrature, accuracy about 1e-8).
// Same as R ptukey(q, nmeans = k, df, nranges = 1) and scipy.stats.studentized_range.cdf(q, k, df).
// qtukey is found by bracketing + bisection/secant on ptukey (R uses a secant iteration; both reach ~1e-9).
import { normCdf, lnGamma } from './stats.js';

const XLEG = [0.981560634246719250690549090149, 0.904117256370474856678465866119, 0.769902674194304687036893833213,
  0.587317954286617447296702418941, 0.367831498998180193752691536644, 0.125233408511468915472441369464];
const ALEG = [0.047175336386511827194615961485, 0.106939325995318430960254718194, 0.160078328543346226334652529543,
  0.203167426723065921749064455810, 0.233492536538354808760849898925, 0.249147045813402785000562436043];
const XLEGQ = [0.989400934991649932596154173450, 0.944575023073232576077988415535, 0.865631202387831743880467897712,
  0.755404408355003033895101194847, 0.617876244402643748446671764049, 0.458016777657227386342419442984,
  0.281603550779258913230460501460, 0.950125098376374401853193354250e-1];
const ALEGQ = [0.271524594117540948517805724560e-1, 0.622535239386478928628438369944e-1, 0.951585116824927848099251076022e-1,
  0.124628971255533872052476282192, 0.149595988816576732081501730547, 0.169156519395002538189312079030,
  0.182603415044923588866763667969, 0.189450610455068496285396723208];

const SQRT_2PI = Math.sqrt(2 * Math.PI);
const pnorm = (x, mu = 0) => normCdf(x - mu);

/** Probability integral of Hartley's form of the range: distribution of the range of `cc` normals (rr ranges). */
function wprob(w, rr, cc) {
  const nleg = 12; const ihalf = 6;
  const C1 = -30; const C2 = -50; const C3 = 60; const bb = 8; const wlar = 3; const wincr1 = 2; const wincr2 = 3;
  const qsqz = w * 0.5;
  if (qsqz >= bb) return 1;
  let prW = 2 * pnorm(qsqz) - 1;
  prW = prW >= Math.exp(C2 / cc) ? prW ** cc : 0;
  const wincr = w > wlar ? wincr1 : wincr2;
  let blb = qsqz;
  const binc = (bb - qsqz) / wincr;
  let bub = blb + binc;
  let einsum = 0;
  const cc1 = cc - 1;
  for (let wi = 1; wi <= wincr; wi++) {
    let elsum = 0;
    const a = 0.5 * (bub + blb);
    const b = 0.5 * (bub - blb);
    for (let jj = 1; jj <= nleg; jj++) {
      let j; let xx;
      if (ihalf < jj) { j = nleg - jj + 1; xx = XLEG[j - 1]; } else { j = jj; xx = -XLEG[j - 1]; }
      const c = b * xx;
      const ac = a + c;
      const qexpo = ac * ac;
      if (qexpo > C3) break;
      const pplus = 2 * pnorm(ac);
      const pminus = 2 * pnorm(ac, w);
      let rinsum = pplus * 0.5 - pminus * 0.5;
      if (rinsum >= Math.exp(C1 / cc1)) {
        rinsum = ALEG[j - 1] * Math.exp(-(0.5 * qexpo)) * rinsum ** cc1;
        elsum += rinsum;
      }
    }
    elsum *= ((2 * b) * cc) / SQRT_2PI;
    einsum += elsum;
    blb = bub;
    bub += binc;
  }
  prW += einsum;
  if (prW <= Math.exp(C1 / rr)) return 0;
  prW **= rr;
  return prW >= 1 ? 1 : prW;
}

/** P(Q ≤ q) for the studentized range of `k` means with `df` error degrees of freedom (R: ptukey(q, k, df)). */
export function ptukey(q, k, df) {
  const nlegq = 16; const ihalfq = 8;
  const eps1 = -30; const eps2 = 1e-14; const dhaf = 100; const dquar = 800; const deigh = 5000; const dlarg = 25000;
  if (Number.isNaN(q) || Number.isNaN(k) || Number.isNaN(df)) return NaN;
  if (q <= 0) return 0;
  if (!(df >= 2) || !(k >= 2)) return NaN;
  if (q === Infinity) return 1;
  if (df > dlarg) return wprob(q, 1, k);
  const f2 = df * 0.5;
  let f2lf = f2 * Math.log(df) - df * Math.LN2 - lnGamma(f2);
  const f21 = f2 - 1;
  const ff4 = df * 0.25;
  let ulen;
  if (df <= dhaf) ulen = 1; else if (df <= dquar) ulen = 0.5; else if (df <= deigh) ulen = 0.25; else ulen = 0.125;
  f2lf += Math.log(ulen);
  let ans = 0;
  for (let i = 1; i <= 50; i++) {
    let otsum = 0;
    const twa1 = (2 * i - 1) * ulen;
    for (let jj = 1; jj <= nlegq; jj++) {
      let j; let t1;
      if (ihalfq < jj) {
        j = jj - ihalfq - 1;
        t1 = f2lf + f21 * Math.log(twa1 + XLEGQ[j] * ulen) - (XLEGQ[j] * ulen + twa1) * ff4;
      } else {
        j = jj - 1;
        t1 = f2lf + f21 * Math.log(twa1 - XLEGQ[j] * ulen) + (XLEGQ[j] * ulen - twa1) * ff4;
      }
      if (t1 >= eps1) {
        const qsqz = ihalfq < jj ? q * Math.sqrt((XLEGQ[j] * ulen + twa1) * 0.5) : q * Math.sqrt((-(XLEGQ[j] * ulen) + twa1) * 0.5);
        otsum += wprob(qsqz, 1, k) * ALEGQ[j] * Math.exp(t1);
      }
    }
    if (i * ulen >= 1 && otsum <= eps2) break;
    ans += otsum;
  }
  return ans > 1 ? 1 : ans;
}

/** Upper tail P(Q ≥ q): the Tukey HSD adjusted p-value. */
export const ptukeySf = (q, k, df) => Math.max(0, 1 - ptukey(q, k, df));

/** Quantile: q with P(Q ≤ q) = p (R: qtukey(p, k, df)); the critical value for Tukey intervals is qtukey(1 − α, k, df). */
export function qtukey(p, k, df) {
  if (Number.isNaN(p) || !(p >= 0 && p <= 1) || !(df >= 2) || !(k >= 2)) return NaN;
  if (p === 0) return 0;
  if (p === 1) return Infinity;
  let lo = 0; let hi = 4;
  while (ptukey(hi, k, df) < p) { lo = hi; hi *= 2; if (hi > 1e4) return NaN; }
  for (let i = 0; i < 200 && hi - lo > 1e-11 * Math.max(1, hi); i++) {
    const mid = 0.5 * (lo + hi);
    if (ptukey(mid, k, df) < p) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}
