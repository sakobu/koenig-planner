import { describe, it, expect } from "vitest";
import { eciToView, rtnToView, type V3 } from "./vec";

// Collapse -0 → 0 so signed zero (rendering-irrelevant) doesn't trip toEqual.
const z = (v: V3): V3 => [v[0] + 0, v[1] + 0, v[2] + 0];
const hypot = (v: V3) => Math.hypot(v[0], v[1], v[2]);
const cross = (a: V3, b: V3): V3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

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

describe("rtnToView", () => {
  it("maps radial (R) to view up (Y)", () => {
    expect(z(rtnToView([1, 0, 0]))).toEqual([0, 1, 0]);
  });

  it("maps transverse (T) to view right (X)", () => {
    expect(z(rtnToView([0, 1, 0]))).toEqual([1, 0, 0]);
  });

  it("maps normal (N) into the screen (view -Z)", () => {
    expect(z(rtnToView([0, 0, 1]))).toEqual([0, 0, -1]);
  });

  it("is the exact basis permutation [t, r, -n]", () => {
    expect(z(rtnToView([1, 2, 3]))).toEqual([2, 1, -3]);
  });

  it("preserves length (proper rotation, no scaling or mirroring)", () => {
    const v: V3 = [0.3, -1.7, 2.4];
    expect(hypot(rtnToView(v))).toBeCloseTo(hypot(v), 12);
  });

  it("is a proper rotation, not a mirror: R×T ↦ N (preserves the deputy's sense of motion)", () => {
    // RTN is right-handed (R×T = N). A proper rotation M obeys M(a)×M(b) = M(a×b);
    // a det -1 mirror would flip the sign, reversing the drawn orbit's direction.
    const R: V3 = [1, 0, 0];
    const T: V3 = [0, 1, 0];
    const N: V3 = [0, 0, 1];
    expect(z(cross(rtnToView(R), rtnToView(T)))).toEqual(z(rtnToView(N)));
  });
});
