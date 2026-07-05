import { useRef, useState } from "react";
import type { ApiError, SolveOutcome, SolveRequest, SolveResponse } from "./wasm";
import { pickDisplay } from "./outcomeDisplay";
import { ErrorBanner } from "./ErrorBanner";
import { Kpis } from "./charts/Kpis";
import { Timeline } from "./charts/Timeline";
import { PrimerMagnitude } from "./charts/PrimerMagnitude";
import { RtnComponents } from "./charts/RtnComponents";
import { PrimerComponents } from "./charts/PrimerComponents";
import { Panel } from "./Panel";
import { PlanTable } from "./PlanTable";
import { EciScene } from "./scene/EciScene";
import { RtnScene } from "./scene/RtnScene";
import { Playback } from "./scene/Playback";

function OkReadout({
  r,
  req,
  error,
}: {
  r: SolveResponse;
  req: SolveRequest;
  error: ApiError | null;
}) {
  const sampleCount = r.geometry.chief_track_eci.length;
  const [index, setIndex] = useState(0);
  // The single clamp for the playback grid. The ChiefGeometry contract keeps all
  // playback-grid arrays (chief_track_eci, deputy_track_rtn, primer_*) equal
  // length, so one clamped frame drives both scenes consistently.
  const frame = Math.min(index, Math.max(0, sampleCount - 1));
  return (
    <section id="output">
      {error && <ErrorBanner kind={error.kind} message={error.message} variant="overlay" />}
      <Kpis r={r} />
      <Panel title="Orbit (ECI)">
        <EciScene g={r.geometry} sampleIndex={frame} />
      </Panel>
      <Panel title="Relative orbit (RTN, chief at origin)">
        <RtnScene g={r.geometry} sampleIndex={frame} />
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
      <Panel title="Plan (precise)">
        <PlanTable req={req} r={r} />
      </Panel>
    </section>
  );
}

export function Readout({
  outcome,
  req,
}: {
  outcome: SolveOutcome | null;
  req: SolveRequest;
}) {
  // Retain the last good response AND the request that produced it: the pairing
  // lets a transient error overlay without unmounting the scenes (keeping camera
  // + scrub state) while the JSON export always ships a response beside the exact
  // request that generated it.
  const lastGood = useRef<{ req: SolveRequest; r: SolveResponse } | null>(null);
  if (outcome?.status === "ok") lastGood.current = { req, r: outcome.value };

  const d = pickDisplay(outcome, lastGood.current?.r ?? null);
  if (d.view === "empty") return <section id="output" />;
  if (d.view === "error") {
    return (
      <section id="output">
        <ErrorBanner kind={d.error.kind} message={d.error.message} />
      </section>
    );
  }
  // The on-screen response is either the current solve (pairs with the current
  // req) or the retained last-good one (pairs with its stored req).
  const shownReq = outcome?.status === "ok" ? req : lastGood.current!.req;
  return <OkReadout r={d.r} req={shownReq} error={d.error} />;
}
