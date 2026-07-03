/** Round a raw step up to a 1/2/5 ×10ⁿ "nice" increment for axis ticks. */
export function niceStep(raw: number): number {
  const exp = Math.floor(Math.log10(raw));
  const f = raw / 10 ** exp;
  const nf = f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10;
  return nf * 10 ** exp;
}

/**
 * Assign a stacking row to each x-position so labels closer than `minGap`
 * (viewBox units) stack onto higher rows instead of overlapping. Processes
 * left-to-right and greedily places each label on the lowest row whose last
 * label is far enough away. Returns the row per input in ORIGINAL order
 * (row 0 = baseline). Lets callers lift colliding labels by `row * lineHeight`.
 */
export function stackRows(xs: number[], minGap: number): number[] {
  const order = xs.map((_, i) => i).sort((a, b) => xs[a] - xs[b]);
  const rows = new Array<number>(xs.length).fill(0);
  const lastX: number[] = []; // last x placed on each row
  for (const i of order) {
    let row = 0;
    while (row < lastX.length && xs[i] - lastX[row] < minGap) row++;
    rows[i] = row;
    lastX[row] = xs[i];
  }
  return rows;
}

/**
 * Largest absolute value over `values`, floored to `floor` so a derived scale
 * stays positive. Spread-free by design: `Math.max(...arr)` on a grid-sized
 * array (>~125k entries) overflows the V8 argument-stack and white-screens.
 */
export function maxAbs(values: number[], floor: number): number {
  let m = floor;
  for (const v of values) {
    const a = Math.abs(v);
    if (a > m) m = a;
  }
  return m;
}
