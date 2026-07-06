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
