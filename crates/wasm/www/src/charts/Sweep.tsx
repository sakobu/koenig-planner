// Cost-vs-horizon trade study: re-solve the plan across a range of final times
// t_f and plot the minimized Δv cost c*(t_f). Longer horizons are generally
// cheaper; the steps are where the optimal burn count changes. The cursor marks
// the current t_f (dragging t_f moves the cursor, not the curve); infeasible
// horizons (solver Err) break the line into gaps.
import { memo } from "react";
import type { SolveRequest } from "../wasm";
import { axisTicks, linePath, niceStep } from "./svgUtil";
import { periodGridTimes } from "../orbit";
import { useSweep, stampWithoutTf, type SweepPoint } from "../useSweep";

const W = 760,
  H = 300;
const padL = 64,
  padR = 30,
  padT = 30,
  padB = 46;

type Feasible = { t_f: number; total_dv: number; nManeuvers: number };
function isFeasible(p: SweepPoint): p is Feasible {
  return "total_dv" in p;
}

export const Sweep = memo(function Sweep({ req, period }: { req: SolveRequest; period: number }) {
  // OkReadout renders only after a successful solve, so wasm is initialized.
  const sweep = useSweep(req, true);
  const stale = sweep.stampHash !== null && sweep.stampHash !== stampWithoutTf(req);
  const feas = sweep.points.filter(isFeasible);

  return (
    <div className={`sweep${stale ? " stale" : ""}`}>
      <div className="sweep-bar">
        <button type="button" onClick={sweep.run} disabled={sweep.status === "running"}>
          {sweep.status === "running"
            ? `Running… ${sweep.done}/${sweep.total}`
            : sweep.stampHash === null
              ? "Run trade study"
              : "Re-run trade study"}
        </button>
        {stale && <span className="sweep-stale">stale — re-run</span>}
      </div>
      {feas.length === 0 ? (
        <p className="sweep-hint">
          {sweep.status === "running"
            ? `solving… ${sweep.done}/${sweep.total}`
            : sweep.stampHash === null
              ? "Re-solve the plan across a range of final times t_f to see how the Δv cost trades against the horizon."
              : "no feasible horizon in range"}
        </p>
      ) : (
        <SweepChart points={sweep.points} feas={feas} tfCurrent={req.t_f} period={period} />
      )}
    </div>
  );
});

function SweepChart({
  points,
  feas,
  tfCurrent,
  period,
}: {
  points: SweepPoint[];
  feas: Feasible[];
  tfCurrent: number;
  period: number;
}) {
  const yBase = H - padB;
  const plotH = yBase - padT;
  const plotW = W - padL - padR;
  const t0 = Math.min(...feas.map((p) => p.t_f));
  const t1 = Math.max(...feas.map((p) => p.t_f));
  const span = Math.max(1e-9, t1 - t0);
  const maxC = feas.reduce((a, p) => Math.max(a, p.total_dv), 0);
  const step = niceStep(Math.max(maxC, 1e-12) / 4);
  const domainMax = Math.max(step, Math.ceil(maxC / step) * step);
  const x = (t: number) => padL + ((t - t0) / span) * plotW;
  const y = (c: number) => yBase - (c / domainMax) * plotH;
  const cursorT = Math.min(Math.max(tfCurrent, t0), t1);

  const yTicks: number[] = [];
  for (let v = 0; v <= domainMax + step / 2; v += step) yTicks.push(v);

  // Contiguous feasible runs → separate polylines (infeasible = gap).
  const runs: Feasible[][] = [];
  let cur: Feasible[] = [];
  for (const p of points) {
    if (isFeasible(p)) cur.push(p);
    else if (cur.length) {
      runs.push(cur);
      cur = [];
    }
  }
  if (cur.length) runs.push(cur);

  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" preserveAspectRatio="xMidYMid meet" className="chart chart-sweep">
      {yTicks.map((v) => (
        <g key={`y${v}`}>
          <line x1={padL} y1={y(v)} x2={W - padR} y2={y(v)} className={v === 0 ? "axis" : "grid"} />
          <text x={padL - 10} y={y(v) + 3.5} className="axis-label" textAnchor="end">{v.toFixed(3)}</text>
        </g>
      ))}
      {axisTicks(t0, t1, 5).map((t) => (
        <text key={`x${t}`} x={x(t)} y={yBase + 18} className="axis-label" textAnchor="middle">{t.toFixed(0)}</text>
      ))}
      {periodGridTimes(t0, t1, period).map((t, k) => (
        <g key={`pg${t}`}>
          <line x1={x(t)} y1={padT} x2={x(t)} y2={yBase} className="period-grid" />
          <text x={x(t)} y={padT - 4} className="mnvr-tag" textAnchor="middle">{`${k + 1}P`}</text>
        </g>
      ))}
      <text x={6} y={15} className="axis-title" textAnchor="start">Δv cost  [m/s]</text>
      <text x={x(t0)} y={yBase + 18} className="axis-label" textAnchor="middle">{t0.toFixed(0)}</text>
      <text x={x(t1)} y={yBase + 18} className="axis-label" textAnchor="middle">{t1.toFixed(0)}</text>
      <text x={padL + plotW / 2} y={yBase + 35} className="axis-title" textAnchor="middle">final time t_f  [s]</text>
      {runs.map((run, i) => (
        <path
          key={i}
          d={linePath(run.map((p) => p.t_f), run.map((p) => p.total_dv), x, y)}
          className="sweep-curve"
        />
      ))}
      {feas.map((p, i) => (
        <circle key={`d${i}`} cx={x(p.t_f)} cy={y(p.total_dv)} r={2.5} className="sweep-dot" />
      ))}
      <line x1={x(cursorT)} y1={padT} x2={x(cursorT)} y2={yBase} className="sweep-cursor" />
      <text x={x(cursorT)} y={padT - 4} className="mnvr-tag" textAnchor="middle">t_f</text>
    </svg>
  );
}
