import { describe, expect, it } from "vitest";
import { clampToWindow, cursorTime } from "./cursorUtil";

describe("cursorTime", () => {
  const times = [0, 30, 60];
  it("returns the scrubbed sample's time", () => {
    expect(cursorTime(times, 1)).toBe(30);
  });
  it("clamps an out-of-range frame into the grid", () => {
    expect(cursorTime(times, 99)).toBe(60);
    expect(cursorTime(times, -1)).toBe(0);
  });
  it("is null without samples", () => {
    expect(cursorTime([], 0)).toBeNull();
  });
  it("snaps to the exact burn time when parked on a burn's nearest sample", () => {
    // grid step 30; burn at 63 rounds to sample index 2 (time 60) — the cursor
    // must report 63 so it lands on the marker drawn at the exact burn time.
    expect(cursorTime([0, 30, 60, 90], 2, [63])).toBe(63);
  });
  it("keeps the grid-sample time away from a burn's sample", () => {
    expect(cursorTime([0, 30, 60, 90], 1, [63])).toBe(30);
  });
});

describe("clampToWindow", () => {
  it("passes a time inside the window", () => {
    expect(clampToWindow(50, 0, 100)).toBe(50);
  });
  it("hides a time outside the window", () => {
    expect(clampToWindow(150, 0, 100)).toBeNull();
    expect(clampToWindow(-1, 0, 100)).toBeNull();
  });
  it("propagates null", () => {
    expect(clampToWindow(null, 0, 100)).toBeNull();
  });
});
