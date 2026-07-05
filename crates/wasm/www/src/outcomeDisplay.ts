import type { ApiError, SolveOutcome, SolveResponse } from "./wasm";

/** What the readout should render for a given outcome. */
export type Display =
  | { view: "empty" }
  | { view: "ok"; r: SolveResponse; error: ApiError | null }
  | { view: "error"; error: ApiError };

/** Decide what the readout shows. On a transient error we keep displaying the
 *  last good response — so the mounted scenes keep their camera and scrub state —
 *  and surface the error as an overlay; only with no prior good solve do we fall
 *  back to a bare error panel. */
export function pickDisplay(
  outcome: SolveOutcome | null,
  lastGood: SolveResponse | null,
): Display {
  if (!outcome) return { view: "empty" };
  if (outcome.status === "ok") return { view: "ok", r: outcome.value, error: null };
  if (lastGood) return { view: "ok", r: lastGood, error: outcome.error };
  return { view: "error", error: outcome.error };
}
