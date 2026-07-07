import { memo } from "react";
import type { SolveResponse } from "../lib/wasm";
import { splitInPlane } from "../lib/export";

export const Kpis = memo(function Kpis({ r }: { r: SolveResponse }) {
  let ipTotal = 0;
  let oopTotal = 0;
  for (const m of r.maneuvers) {
    const { ip, oop } = splitInPlane(m.dv);
    ipTotal += ip;
    oopTotal += oop;
  }
  const cells: [string, string][] = [
    ["Δv cost", `${r.total_dv.toFixed(4)} m/s`],
    ["in-plane Δv", `${ipTotal.toFixed(4)} m/s`],
    ["out-of-plane Δv", `${oopTotal.toFixed(4)} m/s`],
    ["maneuvers", String(r.maneuvers.length)],
    ["iterations", String(r.iterations)],
    ["residual", r.residual.toExponential(2)],
  ];
  return (
    <div className="kpis">
      {cells.map(([label, value]) => (
        <div className="kpi" key={label}>
          <span className="k-label">{label}</span>
          <span className="k-value">{value}</span>
        </div>
      ))}
    </div>
  );
});
