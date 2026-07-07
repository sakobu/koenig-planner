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
      // ~522 km sun-synchronous (i≈97.485° is the exact SSO inclination for this a),
      // near-circular; small bounded relative formation; norm2 (no perigee window).
      chief: { a: 6_900e3, e: 0.001, i: 97.485, raan: 90, argp: 0, mean_anom: 0 },
      t_i: 0,
      t_f: 11_400, // ~2 orbits
      dt: 30,
      w_meters: [0, 200, 100, 0, 0, 100],
      cost: { type: "norm2" },
    },
  },
  {
    id: "leo-coelliptic-hold",
    name: "LEO co-elliptic hold",
    req: {
      // ~500 km, ISS-like inclination. Co-elliptic hold (bounded, δa=0): δe_x
      // traces a coplanar 2:1 R/T ellipse centred ~400 m ahead of the chief —
      // a station-keeping geometry, not a drive-to-docking rendezvous.
      chief: { a: 6_878e3, e: 0.0005, i: 51.6, raan: 0, argp: 0, mean_anom: 0 },
      t_i: 0,
      t_f: 11_400,
      dt: 30,
      w_meters: [0, 400, 150, 0, 0, 0],
      cost: { type: "norm2" },
    },
  },
  {
    id: "geo-relative-slot",
    name: "GEO relative slot offset",
    req: {
      // Geostationary, near-equatorial; small relative longitude/inclination offset
      // (δλ=5 km ≈ 0.007°, ~7% of a 0.1° slot) over half a sidereal day — a
      // station-keeping / formation geometry, not a full slot relocation.
      chief: { a: 42_164e3, e: 0.0002, i: 0.05, raan: 0, argp: 0, mean_anom: 0 },
      t_i: 0,
      t_f: 43_082,
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
