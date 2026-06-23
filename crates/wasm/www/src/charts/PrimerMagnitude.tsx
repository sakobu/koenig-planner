import type { SolveResponse } from "../wasm";

export function PrimerMagnitude({ r }: { r: SolveResponse }) {
  const W = 760,
    H = 300;
  const padL = 58,
    padR = 30,
    padT = 40,
    padB = 44;
  const yBase = H - padB;
  const plotH = yBase - padT;
  const plotW = W - padL - padR;

  const times = r.primer_times;
  const mags = r.primer_magnitude;
  const n = times.length;
  const t0 = n ? times[0] : 0;
  const t1 = n ? times[n - 1] : 1;
  const span = Math.max(1e-9, t1 - t0);
  const peak = mags.reduce((a, b) => Math.max(a, b), 1.0);
  const domainMax = Math.max(1.12, peak + 0.08);
  const x = (t: number) => padL + ((t - t0) / span) * plotW;
  const y = (m: number) => yBase - (m / domainMax) * plotH;

  const path = times.map((t, k) => `${k === 0 ? "M" : "L"}${x(t).toFixed(2)},${y(mags[k]).toFixed(2)}`).join(" ");

  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" preserveAspectRatio="xMidYMid meet" className="chart chart-primer">
      {[0, 0.25, 0.5, 0.75].map((v) => (
        <g key={v}>
          <line x1={padL} y1={y(v)} x2={W - padR} y2={y(v)} className={v === 0 ? "axis" : "grid"} />
          <text x={padL - 10} y={y(v) + 3.5} className="axis-label" textAnchor="end">{v.toFixed(2)}</text>
        </g>
      ))}
      <line x1={padL} y1={y(1)} x2={W - padR} y2={y(1)} className="primer-ref" />
      <text x={W - padR} y={y(1) - 5} className="axis-label" textAnchor="end">|p| = 1</text>
      <text x={6} y={15} className="axis-title" textAnchor="start">primer |p|</text>
      <text x={x(t0)} y={yBase + 18} className="axis-label" textAnchor="middle">{t0.toFixed(0)}</text>
      <text x={x(t1)} y={yBase + 18} className="axis-label" textAnchor="middle">{t1.toFixed(0)}</text>
      <text x={padL + plotW / 2} y={yBase + 35} className="axis-title" textAnchor="middle">time  [s]</text>
      {r.maneuvers.map((m, j) => (
        <line key={j} x1={x(m.t)} y1={padT} x2={x(m.t)} y2={yBase} className="primer-mnvr" />
      ))}
      <path d={path} className="primer-curve" />
      {r.maneuvers.map((m, j) => {
        const idx = times.findIndex((t) => Math.abs(t - m.t) < 1e-6);
        const g = idx >= 0 ? mags[idx] : 1.0;
        return (
          <g key={`d${j}`}>
            <circle cx={x(m.t)} cy={y(g)} r={4} className="stem-dot" />
            <text x={x(m.t)} y={y(g) - 9} className="mnvr-tag" textAnchor="middle">{`mnvr ${j + 1}`}</text>
          </g>
        );
      })}
    </svg>
  );
}
