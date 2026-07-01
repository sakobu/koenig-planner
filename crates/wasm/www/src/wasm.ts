import init, {
  solve,
  version,
  type SolveRequest,
  type SolveResponse,
  type SolveOutcome,
  type CostSpec,
  type ChiefGeometry,
} from "koenig-damico-planner-wasm";

let ready: Promise<void> | null = null;

/** Initialize the wasm module exactly once. */
export function initWasm(): Promise<void> {
  if (!ready) ready = init().then(() => undefined);
  return ready;
}

export { solve, version };
export type {
  SolveRequest,
  SolveResponse,
  SolveOutcome,
  CostSpec,
  ChiefGeometry,
};
