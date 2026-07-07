import { memo } from "react";
import type { SolveResponse } from "../lib/wasm";
import { axisTicks, linePath, nearestIndex, stackRows } from "../lib/svgUtil";
import { periodGridTimes } from "../lib/orbit";
import { cursorTime } from "./cursorUtil";

const W = 760,
  H = 300;
const padL = 58,
  padR = 30,
  padT = 54, // headroom so close-maneuver labels can stack upward
  padB = 44;
const yBase = H - padB;
const plotH = yBase - padT;
const plotW = W - padL - padR;

/** x-scale over the primer grid — shared by the static body and the cursor. */
function xScale(times: number[]): (t: number) => number {
  const n = times.length;
  const t0 = n ? times[0] : 0;
  const t1 = n ? times[n - 1] : 1;
  const span = Math.max(1e-9, t1 - t0);
  return (t) => padL + ((t - t0) / span) * plotW;
}

/* Static layer, memoized on {r, period}: playback ticks (up to 20/s) redraw
 * only the cursor line above, never these grid-sized path strings. */
const Body = memo(function PrimerMagnitudeBody({ r, period }: { r: SolveResponse; period: number }) {
  const times = r.primer_times;
  const mags = r.primer_magnitude;
  const n = times.length;
  const t0 = n ? times[0] : 0;
  const t1 = n ? times[n - 1] : 1;
  const peak = mags.reduce((a, b) => Math.max(a, b), 1.0);
  const domainMax = Math.max(1.12, peak + 0.08);
  const x = xScale(times);
  const y = (m: number) => yBase - (m / domainMax) * plotH;

  const path = linePath(times, mags, x, y);

  const tTicks = axisTicks(t0, t1, 5);
  const pGrid = periodGridTimes(t0, t1, period);

  return (
    <g>
      {[0, 0.25, 0.5, 0.75].map((v) => (
        <g key={v}>
          <line x1={padL} y1={y(v)} x2={W - padR} y2={y(v)} className={v === 0 ? "axis" : "grid"} />
          <text x={padL - 10} y={y(v) + 3.5} className="axis-label" textAnchor="end">{v.toFixed(2)}</text>
        </g>
      ))}
      {tTicks.map((t) => (
        <g key={`tt${t}`}>
          <line x1={x(t)} y1={padT} x2={x(t)} y2={yBase} className="grid" />
          <text x={x(t)} y={yBase + 18} className="axis-label" textAnchor="middle">{t.toFixed(0)}</text>
        </g>
      ))}
      {pGrid.map((t, k) => (
        <g key={`pg${t}`}>
          <line x1={x(t)} y1={padT} x2={x(t)} y2={yBase} className="period-grid" />
          <text x={x(t)} y={padT - 4} className="mnvr-tag" textAnchor="middle">{`${k + 1}P`}</text>
        </g>
      ))}
      <line x1={padL} y1={y(1)} x2={W - padR} y2={y(1)} className="primer-ref" />
      {/* Anchor the reference label at the LEFT of the line: the top-right
          corner is where late maneuvers (at |p|≈1, near t_f) cluster. */}
      <text x={padL + 4} y={y(1) - 5} className="axis-label" textAnchor="start">|p| = 1</text>
      <text x={6} y={15} className="axis-title" textAnchor="start">primer |p|</text>
      <text x={x(t0)} y={yBase + 18} className="axis-label" textAnchor="middle">{t0.toFixed(0)}</text>
      <text x={x(t1)} y={yBase + 18} className="axis-label" textAnchor="middle">{t1.toFixed(0)}</text>
      <text x={padL + plotW / 2} y={yBase + 35} className="axis-title" textAnchor="middle">time  [s]</text>
      {r.maneuvers.map((m, j) => (
        <line key={j} x1={x(m.t)} y1={padT} x2={x(m.t)} y2={yBase} className="primer-mnvr" />
      ))}
      <path d={path} className="primer-curve" />
      {(() => {
        const rows = stackRows(r.maneuvers.map((m) => x(m.t)), 22);
        return r.maneuvers.map((m, j) => {
          const g = mags[nearestIndex(times, m.t)];
          const tagY = Math.max(12, y(g) - 9 - rows[j] * 13);
          return (
            <g key={`d${j}`}>
              <circle cx={x(m.t)} cy={y(g)} r={4} className="stem-dot" />
              <text x={x(m.t)} y={tagY} className="mnvr-tag" textAnchor="middle">{`M${j + 1}`}</text>
            </g>
          );
        });
      })()}
    </g>
  );
});

export const PrimerMagnitude = memo(function PrimerMagnitude({
  r,
  period,
  frame,
}: {
  r: SolveResponse;
  period: number;
  frame: number;
}) {
  const x = xScale(r.primer_times);
  const ct = cursorTime(r.primer_times, frame);
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" preserveAspectRatio="xMidYMid meet" className="chart chart-primer">
      <Body r={r} period={period} />
      {ct !== null && <line x1={x(ct)} y1={padT} x2={x(ct)} y2={yBase} className="time-cursor" />}
    </svg>
  );
});
