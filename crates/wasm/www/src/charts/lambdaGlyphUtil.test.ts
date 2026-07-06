import { describe, it, expect } from "vitest";
import { ROE_LABELS, lambdaBars } from "./lambdaGlyphUtil";

describe("lambdaBars", () => {
  it("labels components in ROE order and normalizes by the largest |λ|", () => {
    const bars = lambdaBars([0.5, -1, 0.25, 0, 0, 0]);
    expect(bars.map((b) => b.label)).toEqual([...ROE_LABELS]);
    expect(bars[0].frac).toBeCloseTo(0.5, 12); // 0.5 / 1
    expect(bars[1].frac).toBeCloseTo(-1, 12); // −1 / 1 (dominant)
    expect(bars[2].frac).toBeCloseTo(0.25, 12);
    expect(bars[1].value).toBe(-1); // raw value preserved for the label
  });
  it("yields all-zero fracs for an all-zero λ (no NaN)", () => {
    const bars = lambdaBars([0, 0, 0, 0, 0, 0]);
    expect(bars.every((b) => b.frac === 0)).toBe(true);
  });
});
