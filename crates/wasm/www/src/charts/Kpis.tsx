import type { SolveResponse } from "../wasm";

export function Kpis({ r }: { r: SolveResponse }) {
  const cells: [string, string][] = [
    ["Δv cost", `${r.total_dv.toFixed(4)} m/s`],
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
}
