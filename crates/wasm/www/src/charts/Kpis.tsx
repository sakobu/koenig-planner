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
  // The Δv cost is the answer, so it's the hero cell; iterations/residual are
  // solver-health telemetry, demoted (dimmer, set apart by a brighter divider)
  // so they don't compete with the physical readouts.
  const cells: { label: string; value: string; cls: string }[] = [
    { label: "Δv cost", value: `${r.total_dv.toFixed(4)} m/s`, cls: "kpi kpi-hero" },
    { label: "in-plane Δv", value: `${ipTotal.toFixed(4)} m/s`, cls: "kpi" },
    { label: "out-of-plane Δv", value: `${oopTotal.toFixed(4)} m/s`, cls: "kpi" },
    { label: "maneuvers", value: String(r.maneuvers.length), cls: "kpi" },
    { label: "iterations", value: String(r.iterations), cls: "kpi kpi-diag" },
    { label: "residual", value: r.residual.toExponential(2), cls: "kpi kpi-diag" },
  ];
  return (
    <div className="kpis">
      {cells.map(({ label, value, cls }) => (
        <div className={cls} key={label}>
          <span className="k-label">{label}</span>
          <span className="k-value">{value}</span>
        </div>
      ))}
    </div>
  );
});
