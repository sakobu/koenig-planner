import { describe, it, expect } from "vitest";
import { MU, chiefPeriod, periodGridTimes } from "./orbit";

describe("chiefPeriod", () => {
  it("uses the core MU (3.986e14), not the CODATA-precise value", () => {
    expect(MU).toBe(3.986e14);
  });
  it("matches the core mean_motion anchor for a = 25 000 km", () => {
    // src/dynamics/orbit.rs: mean_motion(a = 25 000 km) = 1.5971975457e-04 ⇒ T = 2π/n.
    const expected = (2 * Math.PI) / 1.5971975457e-4; // ≈ 39 338.6 s
    expect(chiefPeriod(25_000e3)).toBeCloseTo(expected, 1);
  });
  it("returns NaN for a non-physical axis", () => {
    expect(Number.isNaN(chiefPeriod(-1))).toBe(true);
  });
});

describe("periodGridTimes", () => {
  it("emits interior period multiples within (t0, t1]", () => {
    expect(periodGridTimes(0, 250, 100)).toEqual([100, 200]);
  });
  it("is empty for a sub-period window", () => {
    expect(periodGridTimes(0, 90, 100)).toEqual([]);
  });
  it("is empty for a non-positive or non-finite period", () => {
    expect(periodGridTimes(0, 1000, 0)).toEqual([]);
    expect(periodGridTimes(0, 1000, NaN)).toEqual([]);
  });
});
