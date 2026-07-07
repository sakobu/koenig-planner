import { describe, it, expect } from "vitest";
import { sweepTfValues, stampWithoutTf } from "./useSweep";
import type { SolveRequest } from "../lib/wasm";

const base = { t_i: 0, t_f: 1000, dt: 100 } as unknown as SolveRequest;

describe("sweepTfValues", () => {
  it("snaps every horizon to a whole multiple of dt, clamps ≥ t_i+dt, dedupes", () => {
    const vals = sweepTfValues(base, 48);
    expect(vals.length).toBeGreaterThan(0);
    expect(vals.every((t) => t % 100 === 0)).toBe(true);
    expect(vals.every((t) => t >= 100)).toBe(true);
    expect(new Set(vals).size).toBe(vals.length);
    expect(Math.max(...vals)).toBeLessThanOrEqual(2000); // t_i + 2·span
  });
  it("returns empty for a degenerate window or dt", () => {
    expect(sweepTfValues({ ...base, t_f: 0 } as SolveRequest, 48)).toEqual([]);
    expect(sweepTfValues({ ...base, dt: 0 } as SolveRequest, 48)).toEqual([]);
  });
});

describe("stampWithoutTf", () => {
  it("is identical across two requests differing only in t_f", () => {
    expect(stampWithoutTf(base)).toBe(stampWithoutTf({ ...base, t_f: 5000 } as SolveRequest));
  });
  it("changes when any non-t_f field changes", () => {
    expect(stampWithoutTf(base)).not.toBe(stampWithoutTf({ ...base, dt: 50 } as SolveRequest));
  });
});
