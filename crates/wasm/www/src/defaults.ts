import type { SolveRequest } from "koenig-damico-planner-wasm";

/** The canonical worked-example fixture (angles in degrees, lengths in metres). */
export const GOLDEN: SolveRequest = {
  chief: { a: 25_000e3, e: 0.7, i: 40, raan: 358, argp: 0, mean_anom: 180 },
  t_i: 0,
  t_f: 117_990,
  dt: 30,
  w_metres: [50, 5000, 100, 100, 0, 400],
  cost: { type: "piecewise" },
};
