import { useEffect, useState } from "react";
import { initWasm, version, type SolveRequest } from "./wasm";
import { GOLDEN } from "./defaults";
import { useSolveOutcome } from "./useSolve";
import { Controls } from "./controls/Controls";

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
        <h1>Koenig-D'Amico Maneuver Planner</h1>
        <span className="version">{ready ? `core v${version()}` : "loading…"}</span>
        <span className="status-lamp">{fault ? "● solve fault" : "● plan nominal"}</span>
      </header>
      <main>
        <Controls req={req} setReq={setReq} />
        <section id="output">
          <pre style={{ gridColumn: "1 / -1", fontSize: "0.7rem", color: "#7c8b9a" }}>
            {outcome
              ? outcome.status === "ok"
                ? `ok — ${outcome.value.maneuvers.length} maneuvers, Δv ${outcome.value.total_dv.toFixed(4)} m/s`
                : `${outcome.error.kind}: ${outcome.error.message}`
              : "solving…"}
          </pre>
        </section>
      </main>
    </>
  );
}
