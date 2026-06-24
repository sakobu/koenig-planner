import type { SolveRequest } from "../wasm";
import {
  setChief,
  setWindow,
  setW,
  setCostType,
  setPiecewise,
  type ChiefKey,
} from "./request";
import { NumberField } from "./NumberField";

const CHIEF: { key: ChiefKey; label: string }[] = [
  { key: "a", label: "a [m]" },
  { key: "e", label: "e" },
  { key: "i", label: "i [deg]" },
  { key: "raan", label: "Ω [deg]" },
  { key: "argp", label: "ω [deg]" },
  { key: "mean_anom", label: "M [deg]" },
];
const W_LABELS = ["δa", "δλ", "δe_x", "δe_y", "δi_x", "δi_y"];

// [min, max, step] per field. Ranges chosen for physical sensibility:
// e guarded < 1; angles full-circle; a from a LEO floor to beyond GEO.
const CHIEF_RANGE: Record<ChiefKey, [number, number, number]> = {
  a: [6_778e3, 50_000e3, 1e3],
  e: [0, 0.95, 0.01],
  i: [0, 180, 0.5],
  raan: [0, 360, 1],
  argp: [0, 360, 1],
  mean_anom: [0, 360, 1],
};
const WINDOW_RANGE: Record<"t_i" | "t_f" | "dt", [number, number, number]> = {
  t_i: [0, 200_000, 100],
  t_f: [0, 500_000, 100],
  dt: [1, 600, 1],
};
// w components [m]: along-track (δλ) is the widest; e/i components are tighter.
const W_RANGE: [number, number, number][] = [
  [-10_000, 10_000, 10], // δa
  [-10_000, 10_000, 10], // δλ
  [-2_000, 2_000, 5], // δe_x
  [-2_000, 2_000, 5], // δe_y
  [-2_000, 2_000, 5], // δi_x
  [-2_000, 2_000, 5], // δi_y
];

export function Controls({
  req,
  setReq,
}: {
  req: SolveRequest;
  setReq: (r: SolveRequest) => void;
}) {
  const opt = (v: string) => (v.trim() === "" ? undefined : Number(v));
  const pw = req.cost.type === "piecewise" ? req.cost : null;

  return (
    <form className="controls" onSubmit={(e) => e.preventDefault()}>
      <details className="section" open>
        <summary>Chief orbit</summary>
        {CHIEF.map(({ key, label }) => (
          <NumberField
            key={key}
            label={label}
            value={req.chief[key]}
            onChange={(v) => setReq(setChief(req, key, v))}
            min={CHIEF_RANGE[key][0]}
            max={CHIEF_RANGE[key][1]}
            step={CHIEF_RANGE[key][2]}
          />
        ))}
      </details>

      <details className="section">
        <summary>Window [s]</summary>
        {(["t_i", "t_f", "dt"] as const).map((key) => (
          <NumberField
            key={key}
            label={key}
            value={req[key]}
            onChange={(v) => setReq(setWindow(req, key, v))}
            min={WINDOW_RANGE[key][0]}
            max={WINDOW_RANGE[key][1]}
            step={WINDOW_RANGE[key][2]}
          />
        ))}
      </details>

      <details className="section">
        <summary>Target w [m]</summary>
        {W_LABELS.map((label, idx) => (
          <NumberField
            key={label}
            label={label}
            value={req.w_metres[idx]}
            onChange={(v) => setReq(setW(req, idx, v))}
            min={W_RANGE[idx][0]}
            max={W_RANGE[idx][1]}
            step={W_RANGE[idx][2]}
          />
        ))}
      </details>

      <details className="section">
        <summary>Cost model</summary>
        <label>
          type
          <select
            value={req.cost.type}
            onChange={(e) =>
              setReq(setCostType(req, e.target.value as SolveRequest["cost"]["type"]))
            }
          >
            <option value="piecewise">piecewise</option>
            <option value="norm2">norm2</option>
            <option value="facemax">facemax</option>
          </select>
        </label>
        {pw && (
          <>
            <label>
              period [s]
              <input
                type="number"
                step="any"
                placeholder="auto"
                value={pw.period ?? ""}
                onChange={(e) => setReq(setPiecewise(req, "period", opt(e.target.value)))}
              />
            </label>
            <label>
              t_perigee0 [s]
              <input
                type="number"
                step="any"
                placeholder="auto"
                value={pw.t_perigee0 ?? ""}
                onChange={(e) =>
                  setReq(setPiecewise(req, "t_perigee0", opt(e.target.value)))
                }
              />
            </label>
          </>
        )}
      </details>
    </form>
  );
}
