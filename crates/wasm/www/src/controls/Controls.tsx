import type { SolveRequest } from "../wasm";
import {
  setChief,
  setWindow,
  setW,
  setCostType,
  setPiecewise,
  type ChiefKey,
} from "./request";

const CHIEF: { key: ChiefKey; label: string }[] = [
  { key: "a", label: "a [m]" },
  { key: "e", label: "e" },
  { key: "i", label: "i [deg]" },
  { key: "raan", label: "Ω [deg]" },
  { key: "argp", label: "ω [deg]" },
  { key: "mean_anom", label: "M [deg]" },
];
const W_LABELS = ["δa", "δλ", "δe_x", "δe_y", "δi_x", "δi_y"];

export function Controls({
  req,
  setReq,
}: {
  req: SolveRequest;
  setReq: (r: SolveRequest) => void;
}) {
  const num = (v: string) => Number(v);
  const opt = (v: string) => (v.trim() === "" ? undefined : Number(v));
  const pw = req.cost.type === "piecewise" ? req.cost : null;

  return (
    <form className="controls" onSubmit={(e) => e.preventDefault()}>
      <fieldset>
        <legend>Chief orbit</legend>
        {CHIEF.map(({ key, label }) => (
          <label key={key}>
            {label}
            <input
              type="number"
              step="any"
              value={req.chief[key]}
              onChange={(e) => setReq(setChief(req, key, num(e.target.value)))}
            />
          </label>
        ))}
      </fieldset>

      <fieldset>
        <legend>Window [s]</legend>
        {(["t_i", "t_f", "dt"] as const).map((key) => (
          <label key={key}>
            {key}
            <input
              type="number"
              step="any"
              value={req[key]}
              onChange={(e) => setReq(setWindow(req, key, num(e.target.value)))}
            />
          </label>
        ))}
      </fieldset>

      <fieldset>
        <legend>Target pseudostate w [m]</legend>
        {W_LABELS.map((label, idx) => (
          <label key={label}>
            {label}
            <input
              type="number"
              step="any"
              value={req.w_metres[idx]}
              onChange={(e) => setReq(setW(req, idx, num(e.target.value)))}
            />
          </label>
        ))}
      </fieldset>

      <fieldset>
        <legend>Cost model</legend>
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
      </fieldset>
    </form>
  );
}
