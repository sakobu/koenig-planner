import { useCallback, useEffect, useState } from "react";
import { initWasm, version, type SolveRequest } from "./wasm";
import { GOLDEN } from "./defaults";
import { useSolveOutcome } from "./useSolve";
import { Controls } from "./controls/Controls";
import { Readout } from "./Readout";
import { ErrorBanner } from "./ErrorBanner";

export default function App() {
  const [ready, setReady] = useState(false);
  const [initError, setInitError] = useState<string | null>(null);
  const [req, setReq] = useState<SolveRequest>(GOLDEN);

  // Re-runnable so the init-failure fork can retry (initWasm drops its cached
  // rejection, so a fresh call re-attempts the load).
  const runInit = useCallback(() => {
    setInitError(null);
    initWasm().then(
      () => setReady(true),
      (e: unknown) => setInitError(e instanceof Error ? e.message : String(e)),
    );
  }, []);
  useEffect(() => {
    runInit();
  }, [runInit]);

  const outcome = useSolveOutcome(req, ready);
  const fault = initError !== null || outcome?.status === "err";

  return (
    <>
      <header className={fault ? "fault" : undefined}>
        <h1>Koenig-D'Amico Impulsive Control Solver</h1>
        <span className="version">
          {initError ? `init failed: ${initError}` : ready ? `core v${version()}` : "loading…"}
        </span>
        <span className="status-lamp">{fault ? "● fault" : "● nominal"}</span>
      </header>
      <p className="about">
        Interactive demo of a finite-difference-verified Rust port of the
        Koenig–D'Amico fuel-optimal impulsive control algorithm — minimum-Δv
        maneuver planning for spacecraft relative orbits (ROEs) under J2.{" "}
        <a href="https://github.com/sakobu/koenig-planner" target="_blank" rel="noopener noreferrer">
          GitHub
        </a>
        {" · "}
        <a
          href="https://ieeexplore.ieee.org/document/9209144"
          target="_blank"
          rel="noopener noreferrer"
        >
          Paper (Koenig &amp; D'Amico, IEEE TAC 2020)
        </a>
      </p>
      <main>
        <Controls req={req} setReq={setReq} />
        {initError ? (
          <section id="output">
            <ErrorBanner variant="internal" message={`wasm init failed: ${initError}`} />
            <button type="button" className="retry" onClick={runInit}>
              retry
            </button>
          </section>
        ) : (
          <Readout outcome={outcome} req={req} />
        )}
      </main>
    </>
  );
}
