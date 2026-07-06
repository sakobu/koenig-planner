// ROE phase-plane triptych: the controlled pseudostate δα(t) as three 2D
// projections. Coast polylines end at each burn's EXACT pre-burn point
// (roe_track[k] − roe_jumps[j]); the amber arrow is the exact B·Δv jump.
// Single-series panes: identity is carried by shape (dot = start, ★ = target,
// arrow = burn) plus the panel caption — no legend needed.
import { memo, useMemo } from "react";
import type { SolveResponse } from "../wasm";
import { linePath, niceStep } from "./svgUtil";
import {
  PANES,
  type PaneSpec,
  burnSampleIndex,
  coastSegments,
  jumpArrows,
  paneExtent,
  fmtTick,
} from "./roePlanesUtil";

const W = 320,
  // H = 312 (not 300) makes plotH == plotW == 240 despite the asymmetric
  // top/bottom padding, so the equal-aspect δe/δi panes render a true 1 m = 1 m.
  H = 312;
const padL = 62,
  padR = 18,
  padT = 30,
  padB = 42;
const plotW = W - padL - padR;
const plotH = H - padT - padB;

function JumpArrow({ x1, y1, x2, y2 }: { x1: number; y1: number; x2: number; y2: number }) {
  // A burn can be invisible in one plane (e.g. a pure-N Δv barely moves δe):
  // skip sub-pixel arrows rather than draw a misleading arrowhead blob.
  if (Math.hypot(x2 - x1, y2 - y1) < 2) return null;
  const ang = Math.atan2(y2 - y1, x2 - x1);
  const hx = (d: number) => x2 - 6 * Math.cos(ang + d);
  const hy = (d: number) => y2 - 6 * Math.sin(ang + d);
  return (
    <g className="roe-jump">
      <line x1={x1} y1={y1} x2={x2} y2={y2} />
      <polygon
        points={`${x2},${y2} ${hx(-0.4).toFixed(2)},${hy(-0.4).toFixed(2)} ${hx(0.4).toFixed(2)},${hy(0.4).toFixed(2)}`}
      />
    </g>
  );
}

/** Four-point star marking the target pseudostate w in this plane. */
function TargetStar({ cx, cy }: { cx: number; cy: number }) {
  const r = 6,
    s = 2.3;
  const d =
    `M${cx},${cy - r} L${cx + s},${cy - s} L${cx + r},${cy} L${cx + s},${cy + s} ` +
    `L${cx},${cy + r} L${cx - s},${cy + s} L${cx - r},${cy} L${cx - s},${cy - s} Z`;
  return <path d={d} className="roe-target" />;
}

function Pane({ r, spec, frame }: { r: SolveResponse; spec: PaneSpec; frame: number }) {
  const g = r.geometry;
  // Heavy geometry (O(grid) polylines and their path strings) cached per
  // {r, spec}; each playback tick re-renders only the now-dot.
  const layer = useMemo(() => {
    const burnIdx = r.maneuvers.map((m) => burnSampleIndex(r.primer_times, m.t));
    const segs = coastSegments(g.roe_track, g.roe_jumps, burnIdx, spec.xi, spec.yi);
    const arrows = jumpArrows(g.roe_track, g.roe_jumps, burnIdx, spec.xi, spec.yi);
    const tgt: [number, number] = [g.target_roe[spec.xi], g.target_roe[spec.yi]];
    const ext = paneExtent(segs, tgt, spec.equalAspect);
    const sx = 1.15 * ext.x; // headroom so extremes don't touch the frame
    const sy = 1.15 * ext.y;
    const x = (v: number) => padL + ((v + sx) / (2 * sx)) * plotW;
    const y = (v: number) => padT + plotH - ((v + sy) / (2 * sy)) * plotH;
    // One "nice" tick per half-axis. niceStep rounds up by strictly less than
    // 2.5×, so niceStep(s/3) < 0.834·s — the tick, its grid line, and its label
    // stay inside the frame with margin for any extent (niceStep(s/2) could round
    // up past s and clip against the viewBox).
    const tx = niceStep(sx / 3);
    const ty = niceStep(sy / 3);
    const jsx = (
      <>
        <text x={12} y={16} className="axis-title" textAnchor="start">
          {spec.yLabel}  [m]
        </text>
        <line x1={padL} y1={y(0)} x2={W - padR} y2={y(0)} className="zero-axis" />
        <line x1={x(0)} y1={padT} x2={x(0)} y2={H - padB} className="zero-axis" />
        {[-tx, tx].map((v) => (
          <g key={`x${v}`}>
            <line x1={x(v)} y1={padT} x2={x(v)} y2={H - padB} className="grid" />
            <text x={x(v)} y={H - padB + 16} className="axis-label" textAnchor="middle">
              {fmtTick(v)}
            </text>
          </g>
        ))}
        {[-ty, ty].map((v) => (
          <g key={`y${v}`}>
            <line x1={padL} y1={y(v)} x2={W - padR} y2={y(v)} className="grid" />
            <text x={padL - 8} y={y(v) + 3.5} className="axis-label" textAnchor="end">
              {fmtTick(v)}
            </text>
          </g>
        ))}
        <text x={padL + plotW / 2} y={H - padB + 34} className="axis-title" textAnchor="middle">
          {spec.xLabel}  [m]
        </text>
        {segs.map((pts, i) => (
          <path
            key={i}
            d={linePath(
              pts.map((p) => p[0]),
              pts.map((p) => p[1]),
              x,
              y,
            )}
            className="roe-track"
          />
        ))}
        {arrows.map((a, j) => (
          <JumpArrow key={j} x1={x(a.from[0])} y1={y(a.from[1])} x2={x(a.to[0])} y2={y(a.to[1])} />
        ))}
        <circle cx={x(0)} cy={y(0)} r={4} className="roe-start" />
        <TargetStar cx={x(tgt[0])} cy={y(tgt[1])} />
      </>
    );
    return { x, y, jsx };
  }, [r, g, spec]);
  const now = g.roe_track[Math.min(frame, g.roe_track.length - 1)];
  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="roe-pane" preserveAspectRatio="xMidYMid meet">
      {layer.jsx}
      {now && <circle cx={layer.x(now[spec.xi])} cy={layer.y(now[spec.yi])} r={3.5} className="roe-now" />}
    </svg>
  );
}

export const RoePlanes = memo(function RoePlanes({ r, frame }: { r: SolveResponse; frame: number }) {
  if (r.geometry.roe_track.length === 0) {
    return (
      <div className="chart chart-roeplanes">
        <svg viewBox={`0 0 ${W} ${H}`} className="roe-pane" preserveAspectRatio="xMidYMid meet">
          <text x={W / 2} y={H / 2} className="axis-label" textAnchor="middle">
            trajectory unavailable
          </text>
        </svg>
      </div>
    );
  }
  return (
    <div className="chart chart-roeplanes">
      {PANES.map((spec) => (
        <Pane key={spec.key} r={r} spec={spec} frame={frame} />
      ))}
    </div>
  );
});
