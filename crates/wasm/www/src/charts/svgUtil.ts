/** Build an SVG `path` `d` from parallel `xs`/`ys` sample arrays mapped through
 *  the plot's `x`/`y` scales: `M` to the first point, `L` to the rest. Empty
 *  when there are no samples. Coordinates are rounded to 2 dp to keep the path
 *  string compact. */
export function linePath(
  xs: number[],
  ys: number[],
  x: (v: number) => number,
  y: (v: number) => number,
): string {
  return xs
    .map((xv, k) => `${k === 0 ? "M" : "L"}${x(xv).toFixed(2)},${y(ys[k]).toFixed(2)}`)
    .join(" ");
}

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

/** Index of the sample in `times` nearest to `t` (linear scan; `times` is a
 *  monotonic grid in practice). Robust to a burn time that falls between grid
 *  samples: always returns a real sample index, never an out-of-range sentinel. */
export function nearestIndex(times: number[], t: number): number {
  let best = 0;
  let bestD = Infinity;
  for (let k = 0; k < times.length; k++) {
    const d = Math.abs(times[k] - t);
    if (d < bestD) {
      bestD = d;
      best = k;
    }
  }
  return best;
}

/** "Nice" interior tick values strictly between `t0` and `t1` — about `count`
 *  of them, on a 1/2/5×10ⁿ step. Endpoints are excluded (drawn separately).
 *  Empty when the span is non-positive. */
export function axisTicks(t0: number, t1: number, count: number): number[] {
  if (!(t1 - t0 > 0)) return [];
  const step = niceStep((t1 - t0) / count);
  const out: number[] = [];
  for (let t = Math.ceil(t0 / step) * step; t < t1 - 1e-9; t += step) {
    if (t > t0 + 1e-9) out.push(t);
  }
  return out;
}
