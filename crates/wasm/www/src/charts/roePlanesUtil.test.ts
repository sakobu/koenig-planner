import { describe, it, expect } from "vitest";
import {
  PANES,
  burnSampleIndex,
  coastSegments,
  jumpArrows,
  paneExtent,
  fmtTick,
} from "./roePlanesUtil";

// A 4-sample track with one burn landing on sample 2 (δe pane components 2/3
// jump by 30/40 there; samples 0-1 are pre-burn zeros).
const track = [
  [0, 0, 0, 0, 0, 0],
  [0, 0, 0, 0, 0, 0],
  [10, 20, 30, 40, 50, 60],
  [11, 21, 31, 41, 51, 61],
];
const jumps = [[10, 20, 30, 40, 50, 60]];

describe("PANES", () => {
  it("projects the documented component pairs (δe, δi, δa–δλ)", () => {
    expect(PANES.map((p) => [p.xi, p.yi])).toEqual([
      [2, 3],
      [4, 5],
      [1, 0],
    ]);
    expect(PANES.map((p) => p.equalAspect)).toEqual([true, true, false]);
  });
});

describe("burnSampleIndex", () => {
  it("finds the exact grid sample", () => {
    expect(burnSampleIndex([0, 30, 60, 90], 60)).toBe(2);
  });
  it("falls back to the nearest sample", () => {
    expect(burnSampleIndex([0, 30, 60, 90], 44)).toBe(1);
  });
});

describe("coastSegments", () => {
  it("ends the pre-burn coast at the exact pre-burn point", () => {
    const segs = coastSegments(track, jumps, [2], 2, 3);
    expect(segs).toHaveLength(2);
    // samples 0,1 then the exact pre-burn point track[2] − jumps[0] = (0, 0)
    expect(segs[0]).toEqual([
      [0, 0],
      [0, 0],
      [0, 0],
    ]);
    // post-burn sample opens the next coast
    expect(segs[1]).toEqual([
      [30, 40],
      [31, 41],
    ]);
  });
  it("degenerates gracefully when the burn is at sample 0", () => {
    const t2 = [
      [5, 5, 5, 5, 5, 5],
      [6, 6, 6, 6, 6, 6],
    ];
    const segs = coastSegments(t2, [[5, 5, 5, 5, 5, 5]], [0], 0, 1);
    expect(segs[0]).toEqual([[0, 0]]);
    expect(segs[1]).toEqual([
      [5, 5],
      [6, 6],
    ]);
  });
});

describe("jumpArrows", () => {
  it("runs from the exact pre-burn point to the post-burn sample", () => {
    const [a] = jumpArrows(track, jumps, [2], 2, 3);
    expect(a).toEqual({ from: [0, 0], to: [30, 40] });
  });
});

describe("paneExtent", () => {
  it("covers track and target, symmetric about zero", () => {
    const ext = paneExtent([[[-20, 5], [10, -8]]], [15, 30], false);
    expect(ext).toEqual({ x: 20, y: 30 });
  });
  it("couples both axes when equal-aspect", () => {
    const ext = paneExtent([[[-20, 5]]], [15, 30], true);
    expect(ext).toEqual({ x: 30, y: 30 });
  });
  it("floors the extent so an all-zero pane still has a frame", () => {
    expect(paneExtent([[[0, 0]]], [0, 0], true)).toEqual({ x: 1, y: 1 });
  });
});

describe("fmtTick", () => {
  it("prints whole meters at chart scale and 2 significant digits below 10 m", () => {
    expect(fmtTick(-5000)).toBe("-5000");
    expect(fmtTick(50)).toBe("50");
    expect(fmtTick(0.5)).toBe("0.50");
  });
});
