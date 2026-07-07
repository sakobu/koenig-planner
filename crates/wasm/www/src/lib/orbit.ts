/** Orbital-mechanics helpers for the demo, kept faithful to the core solver. */

/** Earth gravitational parameter [m³/s²]. Deliberately the CORE's value
 *  (`src/dynamics/constants.rs`: `MU = 3.986e14`), NOT the CODATA-precise
 *  3.986004418e14 — so a client-side period matches the solver's own
 *  `mean_motion = √(MU/a³)` (`src/dynamics/orbit.rs`) exactly. */
export const MU = 3.986e14;

/** Chief orbital period [s] from semimajor axis `a` [m]: 2π·√(a³/μ) = TAU / n.
 *  Non-finite / non-positive `a` yields NaN (the callers guard on that). */
export function chiefPeriod(a: number): number {
  return 2 * Math.PI * Math.sqrt((a * a * a) / MU);
}

/** Absolute grid times `t0 + k·period` (k = 1, 2, …) within `(t0, t1]`, for
 *  chief-period gridlines. Empty when `period` is not a positive finite number
 *  or the window spans less than one period. */
export function periodGridTimes(t0: number, t1: number, period: number): number[] {
  if (!Number.isFinite(period) || period <= 0) return [];
  const out: number[] = [];
  for (let t = t0 + period; t <= t1 + 1e-9; t += period) out.push(t);
  return out;
}
