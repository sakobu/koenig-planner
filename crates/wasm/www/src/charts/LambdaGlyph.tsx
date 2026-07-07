// λ dual certificate: the optimal KKT multiplier as six zero-centered signed
// bars. The primer p(t) = Γᵀ(t)·λ, so the magnitude/component panels are this
// vector projected into control space over time. Bars are normalized to the
// largest |component|; the raw value is printed at each tip.
import { memo } from "react";
import type { SolveResponse } from "../lib/wasm";
import { lambdaBars } from "./lambdaGlyphUtil";

export const LambdaGlyph = memo(function LambdaGlyph({ r }: { r: SolveResponse }) {
  const bars = lambdaBars(r.lambda);
  const W = 760,
    rowH = 26,
    padT = 14,
    padB = 14,
    padL = 62,
    padR = 96;
  const H = padT + bars.length * rowH + padB;
  const half = (W - padL - padR) / 2;
  const cx = padL + half; // zero line
  const bh = 12;
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" preserveAspectRatio="xMidYMid meet" className="chart chart-lambda">
      <line x1={cx} y1={padT} x2={cx} y2={H - padB} className="zero-axis" />
      {bars.map((b, i) => {
        const yc = padT + i * rowH + rowH / 2;
        const len = Math.abs(b.frac) * half;
        const pos = b.frac >= 0;
        return (
          <g key={b.label}>
            <text x={14} y={yc + 4} className="row-label" textAnchor="start">
              {b.label}
            </text>
            <rect
              x={pos ? cx : cx - len}
              y={yc - bh / 2}
              width={Math.max(len, 0.75)}
              height={bh}
              rx={1.5}
              className="lambda-bar"
            />
            <text x={W - padR + 8} y={yc + 4} className="val-label" textAnchor="start">
              {`${pos ? "+" : "−"}${Math.abs(b.value).toFixed(4)}`}
            </text>
          </g>
        );
      })}
    </svg>
  );
});
