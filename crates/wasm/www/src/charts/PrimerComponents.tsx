// Primer vector p(t) = Γᵀ(t)·λ, R/T/N components over time.
// This is the primer (the dual mapped into control space), not the executed
// thrust direction — the optimal impulse fires along the support image s(Γᵀλ),
// parallel to the primer only for the L2 cost. Reading it alongside the
// magnitude panel shows which way the dual rewards thrust at each time.
import { memo } from "react";
import type { SolveResponse } from "../wasm";
import { maxAbs } from "./svgUtil";

const RTN_COLORS = { R: "#ff6b6b", T: "#4dd2ff", N: "#ffb454" } as const;
const RTN_NAME = { R: "radial", T: "transverse", N: "normal" } as const;

function pathFor(times: number[], ys: number[], x: (t: number) => number, y: (v: number) => number): string {
  if (times.length === 0) return "";
  return times.map((t, k) => `${k === 0 ? "M" : "L"}${x(t).toFixed(2)},${y(ys[k]).toFixed(2)}`).join(" ");
}

export const PrimerComponents = memo(function PrimerComponents({ r }: { r: SolveResponse }) {
  const W = 760,
    H = 280;
  const padL = 58,
    padR = 30,
    padT = 40,
    padB = 44;

  const times = r.primer_times;
  const rtn = r.primer_rtn;
  const n = times.length;
  const t0 = n ? times[0] : 0;
  const t1 = n ? times[n - 1] : 1;
  const span = Math.max(1e-9, t1 - t0);
  const maxComp = maxAbs(rtn.flat(), 1e-12);
  const domainMax = maxComp * 1.1;
  const cy0 = padT + (H - padT - padB) / 2; // zero axis
  const half = (H - padT - padB) / 2;
  const x = (t: number) => padL + ((t - t0) / span) * (W - padL - padR);
  const y = (v: number) => cy0 - (v / domainMax) * half;

  // Legend x-positions spread across the header row.
  const lstep = (W - padR - 120 - padL) / 2;

  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      width="100%"
      preserveAspectRatio="xMidYMid meet"
      className="chart chart-primer-rtn"
    >
      {/* Legend */}
      {(["R", "T", "N"] as const).map((comp, k) => {
        const lx = padL + k * lstep;
        return (
          <g key={comp}>
            <rect x={lx} y={padT - 30} width={11} height={11} rx={2} fill={RTN_COLORS[comp]} />
            <text x={lx + 17} y={padT - 20} className="legend-label">
              {RTN_NAME[comp]}
            </text>
          </g>
        );
      })}

      {/* Centered zero axis */}
      <line x1={padL} y1={cy0} x2={W - padR} y2={cy0} className="zero-axis" />
      <text x={padL - 10} y={cy0 + 3.5} className="axis-label" textAnchor="end">
        0
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

      {/* Maneuver vertical guides */}
      {r.maneuvers.map((m, j) => (
        <line key={j} x1={x(m.t)} y1={padT} x2={x(m.t)} y2={H - padB} className="primer-mnvr" />
      ))}

      {/* Three component traces (R/T/N) */}
      {(["R", "T", "N"] as const).map((comp, k) => {
        const ys = rtn.map((p) => p[k]);
        return (
          <path
            key={comp}
            d={pathFor(times, ys, x, y)}
            className="primer-comp"
            stroke={RTN_COLORS[comp]}
          />
        );
      })}
    </svg>
  );
});
