/** Round a raw step up to a 1/2/5 ×10ⁿ "nice" increment for axis ticks. */
export function niceStep(raw: number): number {
  const exp = Math.floor(Math.log10(raw));
  const f = raw / 10 ** exp;
  const nf = f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10;
  return nf * 10 ** exp;
}
