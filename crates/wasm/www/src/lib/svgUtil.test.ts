import { describe, it, expect } from "vitest";
import { nearestIndex, axisTicks } from "./svgUtil";

describe("nearestIndex", () => {
  it("finds the exact grid sample", () => {
    expect(nearestIndex([0, 30, 60, 90], 60)).toBe(2);
  });
  it("falls back to the nearest sample when off-grid", () => {
    expect(nearestIndex([0, 30, 60, 90], 44)).toBe(1);
  });
});

describe("axisTicks", () => {
  it("emits interior 1/2/5 ticks strictly between the endpoints", () => {
    expect(axisTicks(0, 100, 5)).toEqual([20, 40, 60, 80]);
  });
  it("is empty when the endpoints are too close for a step", () => {
    expect(axisTicks(0, 0, 5)).toEqual([]);
  });
});
