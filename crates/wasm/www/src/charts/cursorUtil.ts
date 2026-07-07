/** Scrub-cursor helpers shared by the time charts. */
import { nearestIndex } from "../lib/svgUtil";

/** Time of the scrubbed sample, with the frame clamped into the grid; null when
 *  there are no samples. When `burnTimes` is given and the frame is parked on a
 *  burn's nearest sample, report the EXACT burn time instead of the quantized
 *  grid time — so a step-to-burn lands the cursor on the marker (drawn at `m.t`)
 *  rather than up to ±½ a grid step away from it. */
export function cursorTime(times: number[], frame: number, burnTimes?: number[]): number | null {
  if (times.length === 0) return null;
  const i = Math.min(Math.max(0, frame), times.length - 1);
  if (burnTimes) {
    for (const bt of burnTimes) if (nearestIndex(times, bt) === i) return bt;
  }
  return times[i];
}

/** The cursor time when it falls inside [lo, hi], else null — for charts whose
 *  x-domain is narrower than the playback grid (the burn-window timeline). */
export function clampToWindow(t: number | null, lo: number, hi: number): number | null {
  return t !== null && t >= lo && t <= hi ? t : null;
}
