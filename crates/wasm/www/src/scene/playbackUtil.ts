/** Pure helpers for the playback instrument (time readout, burn ticks on the
 *  scrub track, step-to-burn). Times are absolute grid seconds; tick fractions
 *  are 0..1 in SAMPLE-INDEX space, which is the range input's axis (and equals
 *  time-fraction on the uniform grid). */
import { nearestIndex } from "../charts/svgUtil";

/** "t +32.79 h" — elapsed hours since the horizon epoch `t0`. */
export function fmtHours(t: number, t0: number): string {
  return `t +${((t - t0) / 3600).toFixed(2)} h`;
}

/** "orbit 3.2" — elapsed chief periods since `t0`; em-dash when the period is
 *  not a positive finite number. */
export function fmtOrbit(t: number, t0: number, period: number): string {
  if (!Number.isFinite(period) || period <= 0) return "orbit —";
  return `orbit ${((t - t0) / period).toFixed(1)}`;
}

/** "ν 141°" — true anomaly in degrees on [0, 360). */
export function fmtNu(nuRad: number): string {
  const deg = ((nuRad * 180) / Math.PI) % 360;
  return `ν ${(deg < 0 ? deg + 360 : deg).toFixed(0)}°`;
}

/** Burn positions as 0..1 fractions along the scrub range. Empty when the grid
 *  can't host a slider (fewer than 2 samples). */
export function burnTickFractions(times: number[], burnTimes: number[]): number[] {
  if (times.length < 2) return [];
  return burnTimes.map((t) => nearestIndex(times, t) / (times.length - 1));
}

function burnFrames(times: number[], burnTimes: number[]): number[] {
  return burnTimes.map((t) => nearestIndex(times, t)).sort((a, b) => a - b);
}

/** Sample index of the nearest burn strictly after `frame`; null at/after the
 *  last burn. */
export function nextBurnFrame(times: number[], burnTimes: number[], frame: number): number | null {
  for (const k of burnFrames(times, burnTimes)) if (k > frame) return k;
  return null;
}

/** Sample index of the nearest burn strictly before `frame`; null at/before
 *  the first burn. */
export function prevBurnFrame(times: number[], burnTimes: number[], frame: number): number | null {
  let out: number | null = null;
  for (const k of burnFrames(times, burnTimes)) if (k < frame) out = k;
  return out;
}
