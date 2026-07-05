import type { SolveRequest, SolveResponse } from "./wasm";
import { downloadBlob, toBurnCsv, toPlanJson } from "./export";

/** Full-precision companion to the rounded charts: an export bar (JSON / CSV)
 *  over a monospace burn table that reads raw response numbers — no toFixed. The
 *  root `.plan-table` class is the `#output` grid hook (see style.css). */
export function PlanTable({ req, r }: { req: SolveRequest; r: SolveResponse }) {
  const nu = r.geometry.maneuver_nu;
  return (
    <div className="plan-table">
      <div className="export-bar">
        <button
          type="button"
          onClick={() => downloadBlob("koenig-plan.json", "application/json", toPlanJson(req, r))}
        >
          Download plan (JSON)
        </button>
        <button
          type="button"
          onClick={() => downloadBlob("koenig-burns.csv", "text/csv", toBurnCsv(r))}
        >
          Download burns (CSV)
        </button>
      </div>
      <table className="burns">
        <thead>
          <tr>
            <th>#</th>
            <th>t [s]</th>
            <th>Δv_R</th>
            <th>Δv_T</th>
            <th>Δv_N</th>
            <th>|Δv| [m/s]</th>
            <th>ν [rad]</th>
          </tr>
        </thead>
        <tbody>
          {r.maneuvers.length === 0 ? (
            <tr>
              <td className="empty" colSpan={7}>
                no maneuvers
              </td>
            </tr>
          ) : (
            r.maneuvers.map((m, j) => (
              <tr key={j}>
                <td>{j}</td>
                <td>{m.t}</td>
                <td>{m.dv[0]}</td>
                <td>{m.dv[1]}</td>
                <td>{m.dv[2]}</td>
                <td>{Math.hypot(m.dv[0], m.dv[1], m.dv[2])}</td>
                <td>{nu[j]}</td>
              </tr>
            ))
          )}
        </tbody>
      </table>
    </div>
  );
}
