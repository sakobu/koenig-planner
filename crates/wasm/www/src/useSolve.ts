import { useEffect, useMemo, useRef, useState } from "react";
import { solve, type SolveRequest, type SolveOutcome } from "./wasm";

/** How often, at most, to re-solve while a control is being dragged (ms). */
const SOLVE_INTERVAL_MS = 150;

/**
 * Leading + trailing throttle. Emits `value` immediately once `ms` has elapsed
 * since the last emit, otherwise schedules the latest value for the end of the
 * current window. So the first edit responds instantly and a held drag updates
 * live at ~1/`ms` — versus a debounce, which freezes until the input settles.
 */
export function useThrottled<T>(value: T, ms: number): T {
  const [throttled, setThrottled] = useState(value);
  const lastEmit = useRef(0); // timestamp (ms) of the last emit
  const timer = useRef<number | undefined>(undefined);

  useEffect(() => {
    const emit = () => {
      lastEmit.current = Date.now();
      setThrottled(value);
    };
    const elapsed = Date.now() - lastEmit.current;
    if (elapsed >= ms) {
      emit(); // leading edge
    } else {
      clearTimeout(timer.current);
      timer.current = setTimeout(emit, ms - elapsed); // trailing edge
    }
    return () => clearTimeout(timer.current);
  }, [value, ms]);

  return throttled;
}

/** Solve the (throttled) request once wasm is ready. Memoized on the request.
 *  solve() never throws — it returns the tagged union — and is pure, so the
 *  memo body has no side effects (safe under React 19 StrictMode double-invoke). */
export function useSolveOutcome(
  req: SolveRequest,
  ready: boolean,
): SolveOutcome | null {
  const throttledReq = useThrottled(req, SOLVE_INTERVAL_MS);
  return useMemo(() => (ready ? solve(throttledReq) : null), [throttledReq, ready]);
}
