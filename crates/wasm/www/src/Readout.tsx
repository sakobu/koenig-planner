import type { SolveOutcome } from "./wasm";
import { Kpis } from "./charts/Kpis";
import { Timeline } from "./charts/Timeline";
import { PrimerMagnitude } from "./charts/PrimerMagnitude";
import { Panel } from "./charts/Panel";
import { EciScene } from "./scene/EciScene";

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
  const r = outcome.value;
  return (
    <section id="output">
      <Kpis r={r} />
      <Panel title="Orbit (ECI)">
        <EciScene g={r.geometry} />
      </Panel>
      <Panel title="Δv timeline">
        <Timeline r={r} />
      </Panel>
      <Panel title="Primer magnitude vs time">
        <PrimerMagnitude r={r} />
      </Panel>

    </section>
  );
}
