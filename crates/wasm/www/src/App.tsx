import { useEffect, useState } from "react";
import { initWasm, version } from "./wasm";

export default function App() {
  const [ready, setReady] = useState(false);
  useEffect(() => {
    initWasm().then(() => setReady(true));
  }, []);

  return (
    <>
      <header>
        <h1>Koenig-D'Amico Maneuver Planner</h1>
        <span className="version">{ready ? `core v${version()}` : "loading…"}</span>
      </header>
      <main>{ready ? <p>ready</p> : <p>initializing solver…</p>}</main>
    </>
  );
}
