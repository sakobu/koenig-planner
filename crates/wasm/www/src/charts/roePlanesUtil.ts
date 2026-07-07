/** Pure geometry for the ROE phase-plane triptych. Data space is meters (the
 *  `roe_track` / `target_roe` scaling); the component maps to pixels. */

export type Pt = [number, number];

export interface PaneSpec {
  key: string;
  /** Component indices into a 6-vector `[δa, δλ, δe_x, δe_y, δi_x, δi_y]`. */
  xi: number;
  yi: number;
  xLabel: string;
  yLabel: string;
  /** Equal-aspect panes render 1 m = 1 m on both axes (true geometric planes). */
  equalAspect: boolean;
}

/** δe and δi are true geometric planes (equal aspect: the ω̇ rotation of δe
 *  must render as an arc); δa–δλ pairs different physical quantities, so its
 *  axes scale independently (presets span δλ 200–5000 m vs δa 0–50 m). */
export const PANES: PaneSpec[] = [
  { key: "de", xi: 2, yi: 3, xLabel: "δe_x", yLabel: "δe_y", equalAspect: true },
  { key: "di", xi: 4, yi: 5, xLabel: "δi_x", yLabel: "δi_y", equalAspect: true },
  { key: "dadl", xi: 1, yi: 0, xLabel: "δλ", yLabel: "δa", equalAspect: false },
];

// Nearest-sample lookup lives in svgUtil now (shared with PrimerMagnitude);
// re-exported here under its original name so this module's consumers (RoePlanes)
// are unchanged.
export { nearestIndex as burnSampleIndex } from "../lib/svgUtil";

function proj(v: number[], xi: number, yi: number): Pt {
  return [v[xi], v[yi]];
}

/**
 * Split the track into per-coast polylines (data coords) for one pane. Coast
 * j ends at the EXACT pre-burn point `track[k] − jumps[j]` at the next burn's
 * sample k (the sample itself carries the post-burn value and opens coast
 * j+1), so a jump is never smeared across a grid interval.
 */
export function coastSegments(
  track: number[][],
  jumps: number[][],
  burnIdx: number[],
  xi: number,
  yi: number,
): Pt[][] {
  const segs: Pt[][] = [];
  let start = 0;
  for (let j = 0; j < burnIdx.length; j++) {
    const k = burnIdx[j];
    const pts: Pt[] = [];
    for (let s = start; s < k; s++) pts.push(proj(track[s], xi, yi));
    pts.push([track[k][xi] - jumps[j][xi], track[k][yi] - jumps[j][yi]]);
    segs.push(pts);
    start = k;
  }
  const tail: Pt[] = [];
  for (let s = start; s < track.length; s++) tail.push(proj(track[s], xi, yi));
  segs.push(tail);
  return segs;
}

/** Exact burn arrows: pre-burn point → post-burn sample, per maneuver. */
export function jumpArrows(
  track: number[][],
  jumps: number[][],
  burnIdx: number[],
  xi: number,
  yi: number,
): { from: Pt; to: Pt }[] {
  return burnIdx.map((k, j) => ({
    from: [track[k][xi] - jumps[j][xi], track[k][yi] - jumps[j][yi]] as Pt,
    to: [track[k][xi], track[k][yi]] as Pt,
  }));
}

/**
 * Symmetric-about-zero half-extents covering the coast segments (which already
 * include the pre-burn tails), the target star, and a `floor` (meters) so an
 * all-zero pane keeps a drawable frame. Equal-aspect panes couple both axes.
 */
export function paneExtent(
  segs: Pt[][],
  target: Pt,
  equalAspect: boolean,
  floor = 1,
): { x: number; y: number } {
  let x = Math.abs(target[0]);
  let y = Math.abs(target[1]);
  for (const seg of segs) {
    for (const [px, py] of seg) {
      if (Math.abs(px) > x) x = Math.abs(px);
      if (Math.abs(py) > y) y = Math.abs(py);
    }
  }
  x = Math.max(x, floor);
  y = Math.max(y, floor);
  if (equalAspect) {
    const m = Math.max(x, y);
    return { x: m, y: m };
  }
  return { x, y };
}

/** Compact tick label for meter values: whole meters at chart scale, two
 *  significant digits below 10 m. */
export function fmtTick(v: number): string {
  return Math.abs(v) >= 10 ? v.toFixed(0) : v.toPrecision(2);
}
