import { useCallback, useEffect, useRef, useState } from "react";
import { solve, type SolveRequest } from "./wasm";

export type SweepPoint =
  | { t_f: number; total_dv: number; nManeuvers: number }
  | { t_f: number; feasible: false };

export interface SweepState {
  status: "idle" | "running";
  done: number;
  total: number;
  points: SweepPoint[];
  /** Stamp of the request that produced `points` (null before the first run). */
  stampHash: string | null;
}

/** Horizons sampled per trade study. */
const N_SWEEP = 48;
/** Solves per idle batch — keeps each yield short (~11 ms/solve). */
const BATCH = 4;

/** Stable stamp of everything the sweep depends on EXCEPT `t_f` (its free
 *  variable). Dragging `t_f` only moves the cursor, so it must not invalidate a
 *  computed curve. Key order is stable (the control helpers preserve it). */
export function stampWithoutTf(req: SolveRequest): string {
  return JSON.stringify({ ...req, t_f: null });
}

/** The `t_f` horizons to sample: `n` points across [t_i + span/n, t_i + 2·span]
 *  (span = t_f − t_i), each snapped to a whole multiple of `dt` (clean
 *  commensurate grids) and clamped to ≥ t_i + dt, then de-duplicated. Empty for
 *  a degenerate window or dt. */
export function sweepTfValues(req: SolveRequest, n: number): number[] {
  const { t_i, t_f, dt } = req;
  const span = t_f - t_i;
  if (!(span > 0) || !(dt > 0) || n < 2) return [];
  const lo = t_i + span / n;
  const hi = t_i + 2 * span;
  const seen = new Set<number>();
  const out: number[] = [];
  for (let k = 0; k < n; k++) {
    const raw = lo + ((hi - lo) * k) / (n - 1);
    const snapped = t_i + Math.max(1, Math.round((raw - t_i) / dt)) * dt;
    if (!seen.has(snapped)) {
      seen.add(snapped);
      out.push(snapped);
    }
  }
  return out;
}

/** Yield to the browser between batches so a ~0.5 s sweep never blocks input. */
function schedule(fn: () => void): void {
  if (typeof requestIdleCallback === "function") requestIdleCallback(() => fn());
  else setTimeout(fn, 0);
}

/** On-demand cost-vs-horizon sweep. `run()` re-solves the plan across a range
 *  of `t_f`, chunked across idle callbacks. Superseded by a newer `run()` or by
 *  unmount via a monotonic run id. */
export function useSweep(req: SolveRequest, ready: boolean): SweepState & { run: () => void } {
  const [state, setState] = useState<SweepState>({
    status: "idle",
    done: 0,
    total: 0,
    points: [],
    stampHash: null,
  });
  const runId = useRef(0);

  const run = useCallback(() => {
    if (!ready) return;
    const values = sweepTfValues(req, N_SWEEP);
    const stampHash = stampWithoutTf(req);
    const myRun = ++runId.current;
    const points: SweepPoint[] = [];
    setState({ status: "running", done: 0, total: values.length, points: [], stampHash });

    let k = 0;
    const step = () => {
      if (myRun !== runId.current) return; // superseded or unmounted
      const end = Math.min(k + BATCH, values.length);
      for (; k < end; k++) {
        const t_f = values[k];
        const outcome = solve({ ...req, t_f });
        points.push(
          outcome.status === "ok"
            ? { t_f, total_dv: outcome.value.total_dv, nManeuvers: outcome.value.maneuvers.length }
            : { t_f, feasible: false },
        );
      }
      if (k < values.length) {
        setState((s) => ({ ...s, done: k, points: [...points] }));
        schedule(step);
      } else {
        setState({ status: "idle", done: values.length, total: values.length, points, stampHash });
      }
    };
    schedule(step);
  }, [req, ready]);

  // Cancel any in-flight sweep on unmount (bumps the run id the step checks).
  useEffect(() => () => void ++runId.current, []);

  return { ...state, run };
}
