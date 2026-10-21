// Input checks shared by the public functions. Each throws an Error with an English message that names the argument.

export function needNumber(v, name) {
  if (typeof v !== 'number' || !Number.isFinite(v)) throw new Error(`${name} must be a finite number (got ${String(v)}).`);
  return v;
}
export function needPositive(v, name) {
  needNumber(v, name);
  if (!(v > 0)) throw new Error(`${name} must be greater than 0 (got ${v}).`);
  return v;
}
export function needNonNegative(v, name) {
  needNumber(v, name);
  if (v < 0) throw new Error(`${name} cannot be negative (got ${v}).`);
  return v;
}
export function needInteger(v, name, min = 0) {
  needNumber(v, name);
  if (!Number.isInteger(v) || v < min) throw new Error(`${name} must be a whole number of at least ${min} (got ${v}).`);
  return v;
}
export function needProbability(v, name) {
  needNumber(v, name);
  if (v < 0 || v > 1) throw new Error(`${name} must be between 0 and 1 (got ${v}).`);
  return v;
}
export function needOpenProbability(v, name) {
  needNumber(v, name);
  if (!(v > 0 && v < 1)) throw new Error(`${name} must be strictly between 0 and 1 (got ${v}).`);
  return v;
}
export function needTail(v) {
  if (v !== 'two' && v !== 'greater' && v !== 'less') throw new Error(`tail must be 'two', 'greater' or 'less' (got ${String(v)}).`);
  return v;
}
/** Array of finite numbers with at least `min` entries. */
export function needSample(xs, name, min = 2) {
  if (!Array.isArray(xs)) throw new Error(`${name} must be an array of numbers.`);
  if (xs.length < min) throw new Error(`${name} needs at least ${min} values (got ${xs.length}).`);
  for (let i = 0; i < xs.length; i++) {
    if (typeof xs[i] !== 'number' || !Number.isFinite(xs[i])) throw new Error(`${name}[${i}] is not a finite number (got ${String(xs[i])}).`);
  }
  return xs;
}
