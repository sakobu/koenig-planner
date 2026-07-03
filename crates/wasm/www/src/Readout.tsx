import { useRef, useState } from "react";
import type { SolveOutcome, SolveResponse } from "./wasm";
import { pickDisplay, type ApiError } from "./outcomeDisplay";
import { Kpis } from "./charts/Kpis";
import { Timeline } from "./charts/Timeline";
import { PrimerMagnitude } from "./charts/PrimerMagnitude";
import { RtnComponents } from "./charts/RtnComponents";
import { PrimerComponents } from "./charts/PrimerComponents";
import { Panel } from "./charts/Panel";
import { EciScene } from "./scene/EciScene";
import { RtnScene } from "./scene/RtnScene";
import { Playback } from "./scene/Playback";

function OkReadout({ r, error }: { r: SolveResponse; error: ApiError | null }) {
  const sampleCount = r.geometry.chief_track_eci.length;
  const [index, setIndex] = useState(0);
  return (
    <section id="output">
      {error && (
        <div className={`error ${error.kind} overlay`}>
          {`${error.kind}: ${error.message}`}
        </div>
      )}
      <Kpis r={r} />
      <Panel title="Orbit (ECI)">
        <EciScene g={r.geometry} sampleIndex={Math.min(index, Math.max(0, sampleCount - 1))} />
      </Panel>
      <Panel title="Relative orbit (RTN, chief at origin)">
        <RtnScene g={r.geometry} sampleIndex={Math.min(index, Math.max(0, sampleCount - 1))} />
      </Panel>
      <Panel title="Δv timeline">
        <Timeline r={r} />
      </Panel>
      <Panel title="Primer magnitude vs time">
        <PrimerMagnitude r={r} />
      </Panel>
      <Panel title="Δv components (R/T/N)">
        <RtnComponents r={r} />
      </Panel>
      <Panel title="Primer components (R/T/N)">
        <PrimerComponents r={r} />
      </Panel>
      <Panel title="Playback">
        <Playback count={sampleCount} index={index} setIndex={setIndex} />
      </Panel>
    </section>
  );
}

export function Readout({ outcome }: { outcome: SolveOutcome | null }) {
  // Retain the last good response so a transient error (e.g. a mid-edit
  // bad_request) overlays the error without unmounting the scenes — which would
  // reset their camera poses and the scrub index.
  const lastGood = useRef<SolveResponse | null>(null);
  if (outcome?.status === "ok") lastGood.current = outcome.value;

  const d = pickDisplay(outcome, lastGood.current);
  if (d.view === "empty") return <section id="output" />;
  if (d.view === "error") {
    return (
      <section id="output">
        <div className={`error ${d.error.kind}`}>
          {`${d.error.kind}: ${d.error.message}`}
        </div>
      </section>
    );
  }
  return <OkReadout r={d.r} error={d.error} />;
}
