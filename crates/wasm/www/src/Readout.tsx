import { useRef, useState } from "react";
import type { ApiError, SolveOutcome, SolveRequest, SolveResponse } from "./wasm";
import { pickDisplay } from "./outcomeDisplay";
import { chiefPeriod } from "./orbit";
import { ErrorBanner } from "./ErrorBanner";
import { Kpis } from "./charts/Kpis";
import { LambdaGlyph } from "./charts/LambdaGlyph";
import { Timeline } from "./charts/Timeline";
import { PrimerMagnitude } from "./charts/PrimerMagnitude";
import { RtnComponents } from "./charts/RtnComponents";
import { PrimerComponents } from "./charts/PrimerComponents";
import { RoePlanes } from "./charts/RoePlanes";
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
  const period = chiefPeriod(req.chief.a);
  return (
    <section id="output">
      {error && <ErrorBanner kind={error.kind} message={error.message} variant="overlay" />}
      <Kpis r={r} />
      <Panel
        title="λ dual certificate"
        caption="The optimal dual λ (the KKT certificate). The primer p(t) = Γᵀ(t)·λ, so the magnitude and component panels below are this vector projected into control space over time."
      >
        <LambdaGlyph r={r} />
      </Panel>
      <Panel
        title="Orbit (ECI)"
        caption="Chief orbit and burn geometry in the Earth-centered inertial frame. The amber arc (piecewise cost) is the perigee attitude-constraint window."
      >
        <EciScene g={r.geometry} sampleIndex={frame} />
      </Panel>
      <Panel
        title="Target relative orbit (RTN, chief at origin)"
        caption="The deputy's target relative orbit. Burn markers are schematic anchors on that orbit — only the Δv direction is exact; the true transfer trajectory is a later phase."
      >
        <RtnScene g={r.geometry} sampleIndex={frame} />
      </Panel>
      <Panel
        title="ROE phase planes (δe, δi, δa–δλ)"
        caption="The controlled mean-ROE pseudostate δα(t), accumulated from 0 at t_i: coasts follow the J2 STM, amber arrows are the exact B·Δv jump at each burn, ★ marks the target w. δe and δi panes are equal-aspect."
      >
        <RoePlanes r={r} />
      </Panel>
      <Panel title="Δv timeline" caption="Executed Δv magnitude at each maneuver across the horizon.">
        <Timeline r={r} period={period} />
      </Panel>
      <Panel
        title="Primer magnitude vs time"
        caption="|p(t)| reaches the amber |p| = 1 bound exactly at optimal burn times; touching 1 between burns signals slack in the plan."
      >
        <PrimerMagnitude r={r} period={period} />
      </Panel>
      <Panel
        title="Δv components (R/T/N)"
        caption="Executed Δv per burn in the chief RTN frame — radial (R), along-track (T), cross-track (N)."
      >
        <RtnComponents r={r} />
      </Panel>
      <Panel
        title="Primer components (R/T/N)"
        caption="The primer vector p(t) = Γᵀλ in RTN — the dual certificate; each burn's direction is the support direction of p (parallel to p only under the norm2 cost)."
      >
        <PrimerComponents r={r} period={period} />
      </Panel>
      <Panel title="Playback" caption="Scrub the maneuver grid; both 3D scenes track the selected time.">
        <Playback count={sampleCount} index={index} setIndex={setIndex} />
      </Panel>
      <Panel
        title="Plan (full precision)"
        caption="Full-precision burns (m/s) and downloads. The charts round for display; these values and the exports do not."
      >
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
