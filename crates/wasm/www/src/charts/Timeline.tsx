import { memo } from "react";
import type { SolveResponse } from "../wasm";
import { maxAbs, niceStep, stackRows } from "./svgUtil";

export const Timeline = memo(function Timeline({ r }: { r: SolveResponse }) {
  const W = 760,
    H = 300;
  const padL = 58,
    padR = 30,
    padT = 46,
    padB = 44;
  const yBase = H - padB;
  const plotH = yBase - padT;
  const plotW = W - padL - padR;

  const n = r.maneuvers.length;
  if (n === 0) {
    return (
      <svg viewBox={`0 0 ${W} ${H}`} width="100%" preserveAspectRatio="xMidYMid meet" className="chart chart-timeline">
        <text x={W / 2} y={H / 2} className="axis-label" textAnchor="middle">
          no maneuvers
        </text>
      </svg>
    );
  }

  const mags = r.maneuvers.map((m) => Math.hypot(m.dv[0], m.dv[1], m.dv[2]));
  const maxMag = maxAbs(mags, 1e-12);
  const step = niceStep(maxMag / 4);
  const domainMax = Math.max(step, Math.ceil(maxMag / step) * step);

  const tVals = r.maneuvers.map((m) => m.t);
  const t_i = tVals.length ? Math.min(...tVals) : 0;
  const t_f = tVals.length ? Math.max(...tVals) : 1;
  const inset = 0.1 * plotW;
  const span = Math.max(1e-9, t_f - t_i);
  const x = (t: number) => padL + inset + ((t - t_i) / span) * (plotW - 2 * inset);
  const y = (mag: number) => yBase - (mag / domainMax) * plotH;

  const ticks: number[] = [];
  for (let v = 0; v <= domainMax + step / 2; v += step) ticks.push(v);

  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" preserveAspectRatio="xMidYMid meet" className="chart chart-timeline">
      {ticks.map((v) => (
        <g key={v}>
          <line x1={padL} y1={y(v)} x2={W - padR} y2={y(v)} className={v === 0 ? "axis" : "grid"} />
          <text x={padL - 10} y={y(v) + 3.5} className="axis-label" textAnchor="end">
            {v.toFixed(4)}
          </text>
        </g>
      ))}
      <text x={6} y={15} className="axis-title" textAnchor="start">|Δv|  [m/s]</text>
      <text x={x(t_i)} y={yBase + 18} className="axis-label" textAnchor="middle">{t_i.toFixed(0)}</text>
      <text x={x(t_f)} y={yBase + 18} className="axis-label" textAnchor="middle">{t_f.toFixed(0)}</text>
      <text x={padL + plotW / 2} y={yBase + 35} className="axis-title" textAnchor="middle">burn time  [s]</text>
      {(() => {
        const rows = stackRows(r.maneuvers.map((m) => x(m.t)), 22);
        return r.maneuvers.map((m, j) => {
          const mx = x(m.t);
          const my = y(mags[j]);
          const tagY = Math.max(12, my - 25 - rows[j] * 12);
          return (
            <g key={j}>
              <line x1={mx} y1={yBase} x2={mx} y2={my} className="stem" />
              <circle cx={mx} cy={my} r={4} className="stem-dot" />
              <text x={mx} y={my - 11} className="stem-label" textAnchor="middle">{mags[j].toFixed(4)}</text>
              <text x={mx} y={tagY} className="mnvr-tag" textAnchor="middle">{`M${j + 1}`}</text>
            </g>
          );
        });
      })()}
    </svg>
  );
});
