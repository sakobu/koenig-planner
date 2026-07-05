import type { SolveRequest, SolveResponse } from "./wasm";

/** Full-precision CSV of the burn schedule: one row per maneuver — t_s, the
 *  chief-RTN Δv components (m/s), |Δv| (m/s), and the chief true anomaly at the
 *  burn (rad). Numbers are emitted raw for exact round-tripping; the charts
 *  round for display, the export must not. */
export function toBurnCsv(r: SolveResponse): string {
  const header = "t_s,dv_R,dv_T,dv_N,dv_mag,nu_rad";
  const nu = r.geometry.maneuver_nu;
  const rows = r.maneuvers.map((m, j) => {
    const [dr, dt, dn] = m.dv;
    return `${m.t},${dr},${dt},${dn},${Math.hypot(dr, dt, dn)},${nu[j] ?? ""}`;
  });
  return [header, ...rows].join("\n");
}

/** The whole plan as one reproducible document: the exact request that produced
 *  it plus the complete response (maneuvers, KPIs, primer history, geometry).
 *  Pretty-printed; JSON preserves every f64 exactly. */
export function toPlanJson(req: SolveRequest, r: SolveResponse): string {
  return JSON.stringify({ request: req, response: r }, null, 2);
}

/** Trigger a client-side download of `text` as `filename`: a Blob behind a
 *  transient object URL clicked through a throwaway anchor, with the URL
 *  revoked on a deferred tick so engines that read the blob asynchronously
 *  don't race the revoke. DOM side-effect only, so it is not unit-tested — the
 *  pure builders above are. */
export function downloadBlob(filename: string, mime: string, text: string): void {
  const url = URL.createObjectURL(new Blob([text], { type: mime }));
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}
