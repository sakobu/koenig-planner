import { describe, it, expect } from "vitest";
import { eciToView, type V3 } from "./vec";

// Collapse -0 → 0 so signed zero (rendering-irrelevant) doesn't trip toEqual.
const z = (v: V3): V3 => [v[0] + 0, v[1] + 0, v[2] + 0];
const hypot = (v: V3) => Math.hypot(v[0], v[1], v[2]);

describe("eciToView", () => {
  it("lifts the ECI pole (Z) to view up (Y)", () => {
    expect(z(eciToView([0, 0, 1]))).toEqual([0, 1, 0]);
  });

  it("lays an equatorial point (z=0) flat in the horizontal plane (view y=0)", () => {
    expect(eciToView([3, 4, 0])[1]).toBe(0);
  });

  it("is the exact basis permutation [x, z, -y]", () => {
    expect(z(eciToView([1, 2, 3]))).toEqual([1, 3, -2]);
  });

  it("preserves length (proper rotation, no scaling or mirroring)", () => {
    const v: V3 = [0.3, -1.7, 2.4];
    expect(hypot(eciToView(v))).toBeCloseTo(hypot(v), 12);
  });
});
