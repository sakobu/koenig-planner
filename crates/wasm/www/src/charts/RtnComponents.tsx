import { memo } from "react";
import type { SolveResponse } from "../wasm";
import { niceStep } from "./svgUtil";

const RTN_COLORS = { R: "#ff6b6b", T: "#4dd2ff", N: "#ffb454" } as const;
const RTN_NAME = { R: "radial", T: "transverse", N: "normal" } as const;

export const RtnComponents = memo(function RtnComponents({ r }: { r: SolveResponse }) {
  const n = r.maneuvers.length;

  // Geometry constants — faithful port of the old rtnBars() layout.
  const W = 760,
    rowH = 54,
    padL = 84,
    padR = 104,
    padT = 42,
    padB = 26;
  const H = padT + Math.max(n, 1) * rowH + padB;

  // Guard: if no maneuvers, render a placeholder row.
  if (n === 0) {
    return (
      <svg
        viewBox={`0 0 ${W} ${padT + rowH + padB}`}
        width="100%"
        preserveAspectRatio="xMidYMid meet"
        className="chart chart-rtncomp"
      >
        <text x={W / 2} y={(padT + rowH + padB) / 2 + 4} className="axis-label" textAnchor="middle">
          no maneuvers
        </text>
      </svg>
    );
  }

  const allAbsVals = r.maneuvers.flatMap((m) => m.dv.map(Math.abs));
  const maxComp = Math.max(1e-12, ...allAbsVals);
  const plotW = W - padL - padR;
  const cx = padL + plotW / 2; // zero axis x-coordinate
  const labelRoom = 72; // reserve space at bar tips for value text
  const tickStep = niceStep(maxComp / 3);
  const domainMax = Math.max(tickStep, Math.ceil(maxComp / tickStep) * tickStep);
  const scale = (plotW / 2 - labelRoom) / domainMax;
  const bh = 11,
    gap = 5,
    blockH = 3 * bh + 2 * gap;
  const axisTop = padT - 6;
  const axisBot = padT + n * rowH + 2;

  // Legend x-positions spread across the header row.
  const lstep = (W - padR - 120 - padL) / 2;

  // Vertical grid tick values (symmetric about zero).
  const gridTicks: number[] = [];
  for (let v = -domainMax; v <= domainMax + tickStep / 2; v += tickStep) {
    gridTicks.push(v);
  }

  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      width="100%"
      preserveAspectRatio="xMidYMid meet"
      className="chart chart-rtncomp"
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
      <text x={W - padR} y={padT - 20} className="axis-label" textAnchor="end">
        [m/s]
      </text>

      {/* Vertical grid lines + bottom axis labels */}
      {gridTicks.map((v) => {
        const gx2 = cx + v * scale;
        const isZero = Math.abs(v) < tickStep / 2;
        return (
          <g key={v}>
            <line x1={gx2} y1={axisTop} x2={gx2} y2={axisBot} className={isZero ? "zero-axis" : "grid"} />
            <text x={gx2} y={axisBot + 16} className="axis-label" textAnchor="middle">
              {isZero ? "0" : v.toFixed(3)}
            </text>
          </g>
        );
      })}

      {/* Per-maneuver rows */}
      {r.maneuvers.map((m, j) => {
        const yc = padT + j * rowH + rowH / 2;
        const top = yc - blockH / 2;
        return (
          <g key={j}>
            {/* Row label */}
            <text x={padL - 16} y={yc + 4} className="row-label" textAnchor="end">
              {`mnvr ${j + 1}`}
            </text>
            {/* Three RTN bars */}
            {(["R", "T", "N"] as const).map((comp, k) => {
              const v = m.dv[k];
              const by = top + k * (bh + gap);
              const len = Math.abs(v) * scale;
              const pos = v >= 0;
              return (
                <g key={comp}>
                  <rect
                    x={pos ? cx : cx - len}
                    y={by}
                    width={Math.max(len, 0.75)}
                    height={bh}
                    rx={1.5}
                    fill={RTN_COLORS[comp]}
                  />
                  <text
                    x={pos ? cx + len + 7 : cx - len - 7}
                    y={by + bh - 1.5}
                    className="val-label"
                    textAnchor={pos ? "start" : "end"}
                  >
                    {`${comp} ${pos ? "+" : "−"}${Math.abs(v).toFixed(4)}`}
                  </text>
                </g>
              );
            })}
          </g>
        );
      })}
    </svg>
  );
});
