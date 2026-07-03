import { useEffect, useState } from "react";
import { initWasm, version, type SolveRequest } from "./wasm";
import { GOLDEN } from "./defaults";
import { useSolveOutcome } from "./useSolve";
import { Controls } from "./controls/Controls";
import { Readout } from "./Readout";

export default function App() {
  const [ready, setReady] = useState(false);
  const [initError, setInitError] = useState<string | null>(null);
  const [req, setReq] = useState<SolveRequest>(GOLDEN);
  useEffect(() => {
    initWasm().then(
      () => setReady(true),
      (e: unknown) => setInitError(e instanceof Error ? e.message : String(e)),
    );
  }, []);

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
      <main>
        <Controls req={req} setReq={setReq} />
        {initError ? (
          <section id="output">
            <div className="error internal">{`wasm init failed: ${initError}`}</div>
          </section>
        ) : (
          <Readout outcome={outcome} />
        )}
      </main>
    </>
  );
}
