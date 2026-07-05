import type { SolveRequest } from "../wasm";
import {
  setChief,
  setWindow,
  setW,
  setCostType,
  setPiecewise,
  type ChiefKey,
} from "./request";
import { NumberField, OptionalNumberField } from "./NumberField";

// One descriptor per field: label + [min, max, step], so the label, range, and
// value can't drift out of position. Ranges chosen for physical sensibility:
// e guarded < 1; angles full-circle; a from a LEO floor to beyond GEO.
const CHIEF_FIELDS: { key: ChiefKey; label: string; min: number; max: number; step: number }[] = [
  { key: "a", label: "a [m]", min: 6_778e3, max: 50_000e3, step: 1e3 },
  { key: "e", label: "e", min: 0, max: 0.95, step: 0.01 },
  { key: "i", label: "i [deg]", min: 0, max: 180, step: 0.5 },
  { key: "raan", label: "Ω [deg]", min: 0, max: 360, step: 1 },
  { key: "argp", label: "ω [deg]", min: 0, max: 360, step: 1 },
  { key: "mean_anom", label: "M [deg]", min: 0, max: 360, step: 1 },
];
const WINDOW_FIELDS: { key: "t_i" | "t_f" | "dt"; min: number; max: number; step: number }[] = [
  { key: "t_i", min: 0, max: 200_000, step: 100 },
  { key: "t_f", min: 0, max: 500_000, step: 100 },
  { key: "dt", min: 1, max: 600, step: 1 },
];
// w components [m]: along-track (δλ) is the widest; e/i components are tighter.
const W_FIELDS: { label: string; min: number; max: number; step: number }[] = [
  { label: "δa", min: -10_000, max: 10_000, step: 10 },
  { label: "δλ", min: -10_000, max: 10_000, step: 10 },
  { label: "δe_x", min: -2_000, max: 2_000, step: 5 },
  { label: "δe_y", min: -2_000, max: 2_000, step: 5 },
  { label: "δi_x", min: -2_000, max: 2_000, step: 5 },
  { label: "δi_y", min: -2_000, max: 2_000, step: 5 },
];

export function Controls({
  req,
  setReq,
}: {
  req: SolveRequest;
  setReq: (r: SolveRequest) => void;
}) {
  const pw = req.cost.type === "piecewise" ? req.cost : null;

  return (
    <form className="controls" onSubmit={(e) => e.preventDefault()}>
      <details className="section" open>
        <summary>Chief orbit</summary>
        {CHIEF_FIELDS.map((f) => (
          <NumberField
            key={f.key}
            label={f.label}
            value={req.chief[f.key]}
            onChange={(v) => setReq(setChief(req, f.key, v))}
            min={f.min}
            max={f.max}
            step={f.step}
          />
        ))}
      </details>

      <details className="section">
        <summary>Window [s]</summary>
        {WINDOW_FIELDS.map((f) => (
          <NumberField
            key={f.key}
            label={f.key}
            value={req[f.key]}
            onChange={(v) => setReq(setWindow(req, f.key, v))}
            min={f.min}
            max={f.max}
            step={f.step}
          />
        ))}
      </details>

      <details className="section">
        <summary>Target w [m]</summary>
        {W_FIELDS.map((f, idx) => (
          <NumberField
            key={f.label}
            label={f.label}
            value={req.w_meters[idx]}
            onChange={(v) => setReq(setW(req, idx, v))}
            min={f.min}
            max={f.max}
            step={f.step}
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
            <OptionalNumberField
              label="period [s]"
              value={pw.period}
              onChange={(v) => setReq(setPiecewise(req, "period", v))}
            />
            <OptionalNumberField
              label="t_perigee0 [s]"
              value={pw.t_perigee0}
              onChange={(v) => setReq(setPiecewise(req, "t_perigee0", v))}
            />
          </>
        )}
      </details>
    </form>
  );
}
