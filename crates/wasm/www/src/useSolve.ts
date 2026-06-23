import { useEffect, useMemo, useState } from "react";
import { solve, type SolveRequest, type SolveOutcome } from "./wasm";

/** Debounce any value by `ms` (used to throttle re-solves while dragging). */
export function useDebounced<T>(value: T, ms: number): T {
  const [debounced, setDebounced] = useState(value);
  useEffect(() => {
    const id = setTimeout(() => setDebounced(value), ms);
    return () => clearTimeout(id);
  }, [value, ms]);
  return debounced;
}

/** Solve the (debounced) request once wasm is ready. Memoized on the request.
 *  solve() never throws — it returns the tagged union — and is pure, so the
 *  memo body has no side effects (safe under React 19 StrictMode double-invoke). */
export function useSolveOutcome(
  req: SolveRequest,
  ready: boolean,
): SolveOutcome | null {
  const debouncedReq = useDebounced(req, 150);
  return useMemo(() => (ready ? solve(debouncedReq) : null), [debouncedReq, ready]);
}
