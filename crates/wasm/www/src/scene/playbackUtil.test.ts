import { describe, expect, it } from "vitest";
import {
  burnTickFractions,
  fmtHours,
  fmtNu,
  fmtOrbit,
  nextBurnFrame,
  prevBurnFrame,
} from "./playbackUtil";

const times = [0, 30, 60, 90, 120];

describe("fmtHours", () => {
  it("formats elapsed hours from the epoch", () => {
    expect(fmtHours(118_044, 0)).toBe("t +32.79 h");
    expect(fmtHours(50_000 + 3_600, 50_000)).toBe("t +1.00 h");
  });
});

describe("fmtOrbit", () => {
  it("counts elapsed chief periods", () => {
    expect(fmtOrbit(9_000, 0, 6_000)).toBe("orbit 1.5");
  });
  it("degrades on a non-finite period", () => {
    expect(fmtOrbit(1, 0, NaN)).toBe("orbit —");
  });
});

describe("fmtNu", () => {
  it("renders degrees on [0, 360)", () => {
    expect(fmtNu(Math.PI / 2)).toBe("ν 90°");
    expect(fmtNu(-Math.PI / 2)).toBe("ν 270°");
  });
});

describe("burnTickFractions", () => {
  it("places ticks in sample-index space (the slider's axis)", () => {
    expect(burnTickFractions(times, [60, 120])).toEqual([0.5, 1]);
  });
  it("is empty for a degenerate grid", () => {
    expect(burnTickFractions([0], [0])).toEqual([]);
  });
});

describe("nextBurnFrame / prevBurnFrame", () => {
  it("steps strictly forward to the next burn sample", () => {
    expect(nextBurnFrame(times, [30, 90], 0)).toBe(1);
    expect(nextBurnFrame(times, [30, 90], 1)).toBe(3);
    expect(nextBurnFrame(times, [30, 90], 3)).toBeNull();
  });
  it("steps strictly backward to the previous burn sample", () => {
    expect(prevBurnFrame(times, [30, 90], 4)).toBe(3);
    expect(prevBurnFrame(times, [30, 90], 3)).toBe(1);
    expect(prevBurnFrame(times, [30, 90], 1)).toBeNull();
  });
});
