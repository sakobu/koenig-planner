import { useState } from "react";
import type { SolveOutcome, SolveResponse } from "./wasm";
import { Kpis } from "./charts/Kpis";
import { Timeline } from "./charts/Timeline";
import { PrimerMagnitude } from "./charts/PrimerMagnitude";
import { RtnComponents } from "./charts/RtnComponents";
import { PrimerComponents } from "./charts/PrimerComponents";
import { Panel } from "./charts/Panel";
import { EciScene } from "./scene/EciScene";
import { RtnScene } from "./scene/RtnScene";
import { Playback } from "./scene/Playback";

function OkReadout({ r }: { r: SolveResponse }) {
  const sampleCount = r.geometry.chief_track_eci.length;
  const [index, setIndex] = useState(0);
  return (
    <section id="output">
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
  if (!outcome) return <section id="output" />;
  if (outcome.status === "err") {
    return (
      <section id="output">
        <div className={`error ${outcome.error.kind}`}>
          {`${outcome.error.kind}: ${outcome.error.message}`}
        </div>
      </section>
    );
  }
  return <OkReadout r={outcome.value} />;
}
