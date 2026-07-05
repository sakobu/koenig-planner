import type { SolveRequest } from "./wasm";

/** A named scenario the user can load in one click. Every `req` has been run
 *  through the core solver and returns a non-degenerate plan, and sits inside
 *  the Controls slider ranges. Angles in degrees, lengths in meters.
 *  w_meters order: [δa, δλ, δe_x, δe_y, δi_x, δi_y]. */
export interface Preset {
  id: string;
  name: string;
  req: SolveRequest;
}

export const PRESETS: Preset[] = [
  {
    id: "heo-worked-example",
    name: "HEO worked example",
    req: {
      chief: { a: 25_000e3, e: 0.7, i: 40, raan: 358, argp: 0, mean_anom: 180 },
      t_i: 0,
      t_f: 117_990,
      dt: 30,
      w_meters: [50, 5000, 100, 100, 0, 400],
      cost: { type: "piecewise" },
    },
  },
  {
    id: "leo-sso-formation",
    name: "LEO sun-sync formation",
    req: {
      // ~520 km sun-synchronous, near-circular; norm2 (no perigee window to gauge).
      chief: { a: 6_900e3, e: 0.001, i: 97.4, raan: 90, argp: 0, mean_anom: 0 },
      t_i: 0,
      t_f: 11_400, // ~2 orbits
      dt: 30,
      w_meters: [0, 200, 100, 0, 0, 100],
      cost: { type: "norm2" },
    },
  },
  {
    id: "leo-rendezvous",
    name: "LEO rendezvous (co-elliptic)",
    req: {
      // ~500 km, ISS-like inclination; a pure along-track rephasing.
      chief: { a: 6_878e3, e: 0.0005, i: 51.6, raan: 0, argp: 0, mean_anom: 0 },
      t_i: 0,
      t_f: 11_400,
      dt: 30,
      w_meters: [0, 400, 0, 0, 0, 0],
      cost: { type: "norm2" },
    },
  },
  {
    id: "geo-relocation",
    name: "GEO relocation",
    req: {
      // Geostationary, near-equatorial; a slot change over half a sidereal day.
      chief: { a: 42_164e3, e: 0.0002, i: 0.05, raan: 0, argp: 0, mean_anom: 0 },
      t_i: 0,
      t_f: 43_200,
      dt: 120,
      w_meters: [0, 5000, 500, 0, 0, 500],
      cost: { type: "norm2" },
    },
  },
];

/** The default scenario on first load — the paper's HEO worked example. */
export const GOLDEN: SolveRequest = PRESETS[0].req;

/** The id of the preset whose request exactly equals `req`, or null once the
 *  user has edited away from every preset ("Custom"). Structural compare via
 *  JSON: presets set no optional `params` / `initial_times`, and the control
 *  update helpers preserve key order, so serialized equality is stable here. */
export function presetIdFor(req: SolveRequest): string | null {
  const key = JSON.stringify(req);
  const match = PRESETS.find((p) => JSON.stringify(p.req) === key);
  return match ? match.id : null;
}
