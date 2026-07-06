/** Scrub-cursor helpers shared by the time charts. */

/** Time of the scrubbed sample, with the frame clamped into the grid; null
 *  when there are no samples. */
export function cursorTime(times: number[], frame: number): number | null {
  if (times.length === 0) return null;
  return times[Math.min(Math.max(0, frame), times.length - 1)];
}

/** The cursor time when it falls inside [lo, hi], else null — for charts whose
 *  x-domain is narrower than the playback grid (the burn-window timeline). */
export function clampToWindow(t: number | null, lo: number, hi: number): number | null {
  return t !== null && t >= lo && t <= hi ? t : null;
}
