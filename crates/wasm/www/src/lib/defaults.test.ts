import { describe, it, expect } from "vitest";
import { PRESETS, presetIdFor } from "./defaults";

describe("PRESETS", () => {
  it("has unique ids and non-empty names", () => {
    const ids = PRESETS.map((p) => p.id);
    expect(new Set(ids).size).toBe(ids.length);
    for (const p of PRESETS) expect(p.name.length).toBeGreaterThan(0);
  });

  it("every preset is a structurally valid SolveRequest", () => {
    for (const p of PRESETS) {
      const { chief, t_i, t_f, dt, w_meters, cost } = p.req;
      expect(w_meters).toHaveLength(6);
      for (const w of w_meters) expect(Number.isFinite(w)).toBe(true);
      for (const v of Object.values(chief)) expect(Number.isFinite(v)).toBe(true);
      expect(chief.e).toBeGreaterThanOrEqual(0);
      expect(chief.e).toBeLessThan(1);
      expect(t_f).toBeGreaterThan(t_i);
      expect(dt).toBeGreaterThan(0);
      expect(["norm2", "facemax", "piecewise"]).toContain(cost.type);
    }
  });
});

describe("presetIdFor", () => {
  it("identifies each preset by its exact request", () => {
    for (const p of PRESETS) expect(presetIdFor(p.req)).toBe(p.id);
  });

  it("returns null once a value is edited away from every preset", () => {
    const edited = { ...PRESETS[0].req, dt: PRESETS[0].req.dt + 1 };
    expect(presetIdFor(edited)).toBeNull();
  });
});
