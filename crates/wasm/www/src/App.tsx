import { useEffect, useState } from "react";
import { initWasm, version, type SolveRequest } from "./wasm";
import { GOLDEN } from "./defaults";
import { useSolveOutcome } from "./useSolve";
import { Controls } from "./controls/Controls";
import { Readout } from "./Readout";

export default function App() {
  const [ready, setReady] = useState(false);
  const [req, setReq] = useState<SolveRequest>(GOLDEN);
  useEffect(() => {
    initWasm().then(() => setReady(true));
  }, []);

  const outcome = useSolveOutcome(req, ready);
  const fault = outcome?.status === "err";

  return (
    <>
      <header className={fault ? "fault" : undefined}>
        <h1>Koenig-D'Amico Impulsive Control Solver</h1>
        <span className="version">{ready ? `core v${version()}` : "loading…"}</span>
        <span className="status-lamp">{fault ? "● fault" : "● nominal"}</span>
      </header>
      <main>
        <Controls req={req} setReq={setReq} />
        <Readout outcome={outcome} />
      </main>
    </>
  );
}
