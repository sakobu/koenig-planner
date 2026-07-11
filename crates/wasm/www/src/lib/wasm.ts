import init from "koenig-damico-planner-wasm";

// Curated facade over the generated wasm-bindgen package: the sole importer of
// it, re-exporting only what the app is meant to use — hiding the generated
// noise (raw init/initSync, solve_json, internal DTOs) — behind one local
// specifier that also absorbs package renames.
export { solve, version, sweep_dual } from "koenig-damico-planner-wasm";
export type {
  SolveRequest,
  SolveResponse,
  SolveOutcome,
  SweepRequest,
  SweepOutcome,
  SweepPoint,
  CostSpec,
  ChiefGeometry,
  ApiError,
} from "koenig-damico-planner-wasm";

let ready: Promise<void> | null = null;

/** Initialize the wasm module exactly once. */
export function initWasm(): Promise<void> {
  if (!ready) {
    ready = init()
      .then(() => undefined)
      .catch((e) => {
        // Drop the cached rejection so a later call can retry, rather than
        // permanently returning the failed promise.
        ready = null;
        throw e;
      });
  }
  return ready;
}
