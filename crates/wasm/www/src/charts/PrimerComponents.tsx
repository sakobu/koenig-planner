// Primer vector p(t) = Γᵀ(t)·λ, R/T/N components over time.
// This is the primer (the dual mapped into control space), not the executed
// thrust direction — the optimal impulse fires along the support image s(Γᵀλ),
// parallel to the primer only for the L2 cost. Reading it alongside the
// magnitude panel shows which way the dual rewards thrust at each time.
import { memo } from "react";
import type { SolveResponse } from "../lib/wasm";
import { axisTicks, linePath, maxAbs } from "../lib/svgUtil";
import { periodGridTimes } from "../lib/orbit";
import { RTN_COLORS } from "../lib/rtn";
import { RtnLegend } from "./RtnLegend";
import { cursorTime } from "./cursorUtil";

const W = 760,
  H = 280;
const padL = 58,
  padR = 30,
  padT = 40,
  padB = 44;

/** x-scale over the primer grid — shared by the static body and the cursor. */
function xScale(times: number[]): (t: number) => number {
  const n = times.length;
  const t0 = n ? times[0] : 0;
  const t1 = n ? times[n - 1] : 1;
  const span = Math.max(1e-9, t1 - t0);
  return (t) => padL + ((t - t0) / span) * (W - padL - padR);
}

/* Static layer, memoized on {r, period}: playback ticks (up to 20/s) redraw
 * only the cursor line above, never these grid-sized path strings. */
const Body = memo(function PrimerComponentsBody({ r, period }: { r: SolveResponse; period: number }) {
  const times = r.primer_times;
  const rtn = r.primer_rtn;
  const n = times.length;
  const t0 = n ? times[0] : 0;
  const t1 = n ? times[n - 1] : 1;
  const maxComp = maxAbs(rtn.flat(), 1e-12);
  const domainMax = maxComp * 1.1;
  const cy0 = padT + (H - padT - padB) / 2; // zero axis
  const half = (H - padT - padB) / 2;
  const x = xScale(times);
  const y = (v: number) => cy0 - (v / domainMax) * half;

  const tTicks = axisTicks(t0, t1, 5);
  const pGrid = periodGridTimes(t0, t1, period);

  return (
    <g>
      {/* Legend */}
      <RtnLegend x={padL} y={padT} width={W} padR={padR} />

      {/* Centered zero axis */}
      <line x1={padL} y1={cy0} x2={W - padR} y2={cy0} className="zero-axis" />
      <text x={padL - 10} y={cy0 + 3.5} className="axis-label" textAnchor="end">
        0
      </text>
      <text x={padL - 10} y={y(domainMax) + 3.5} className="axis-label" textAnchor="end">
        {domainMax.toFixed(2)}
      </text>
      <text x={padL - 10} y={y(-domainMax) + 3.5} className="axis-label" textAnchor="end">
        {(-domainMax).toFixed(2)}
      </text>

      {/* Time axis labels */}
      <text x={x(t0)} y={H - padB + 18} className="axis-label" textAnchor="middle">
        {t0.toFixed(0)}
      </text>
      <text x={x(t1)} y={H - padB + 18} className="axis-label" textAnchor="middle">
        {t1.toFixed(0)}
      </text>
      <text x={padL + (W - padL - padR) / 2} y={H - padB + 35} className="axis-title" textAnchor="middle">
        time  [s]
      </text>

      {tTicks.map((t) => (
        <g key={`tt${t}`}>
          <line x1={x(t)} y1={padT} x2={x(t)} y2={H - padB} className="grid" />
          <text x={x(t)} y={H - padB + 18} className="axis-label" textAnchor="middle">{t.toFixed(0)}</text>
        </g>
      ))}
      {pGrid.map((t, k) => (
        <g key={`pg${t}`}>
          <line x1={x(t)} y1={padT} x2={x(t)} y2={H - padB} className="period-grid" />
          <text x={x(t)} y={padT - 4} className="mnvr-tag" textAnchor="middle">{`${k + 1}P`}</text>
        </g>
      ))}

      {/* Maneuver vertical guides, tagged with the shared M1..Mn burn index. */}
      {r.maneuvers.map((m, j) => (
        <g key={j}>
          <line x1={x(m.t)} y1={padT} x2={x(m.t)} y2={H - padB} className="primer-mnvr" />
          <text x={x(m.t)} y={padT + 9} className="mnvr-tag" textAnchor="middle">{`M${j + 1}`}</text>
        </g>
      ))}

      {/* Three component traces (R/T/N) */}
      {(["R", "T", "N"] as const).map((comp, k) => {
        const ys = rtn.map((p) => p[k]);
        return (
          <path
            key={comp}
            d={linePath(times, ys, x, y)}
            className="primer-comp"
            stroke={RTN_COLORS[comp]}
          />
        );
      })}
    </g>
  );
});

export const PrimerComponents = memo(function PrimerComponents({
  r,
  period,
  frame,
}: {
  r: SolveResponse;
  period: number;
  frame: number;
}) {
  const x = xScale(r.primer_times);
  const ct = cursorTime(r.primer_times, frame, r.maneuvers.map((m) => m.t));
  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      width="100%"
      preserveAspectRatio="xMidYMid meet"
      className="chart chart-primer-rtn"
    >
      <Body r={r} period={period} />
      {ct !== null && <line x1={x(ct)} y1={padT} x2={x(ct)} y2={H - padB} className="time-cursor" />}
    </svg>
  );
});
