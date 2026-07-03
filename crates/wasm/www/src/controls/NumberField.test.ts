import { describe, it, expect } from "vitest";
import { parseCommit } from "./NumberField";

describe("parseCommit", () => {
  it("commits a plain number", () => {
    expect(parseCommit("42")).toBe(42);
    expect(parseCommit("6800000")).toBe(6800000);
    expect(parseCommit("0.001")).toBe(0.001);
  });

  it("commits negatives and scientific notation", () => {
    expect(parseCommit("-0.5")).toBe(-0.5);
    expect(parseCommit("6.674e-11")).toBe(6.674e-11);
  });

  it("holds an empty or whitespace draft without committing (no 0-snap)", () => {
    expect(parseCommit("")).toBeNull();
    expect(parseCommit("   ")).toBeNull();
  });

  it("holds an in-progress lone minus without clobbering to NaN", () => {
    expect(parseCommit("-")).toBeNull();
  });

  it("holds an in-progress bare exponent", () => {
    expect(parseCommit("1e")).toBeNull();
  });

  it("rejects non-numeric input", () => {
    expect(parseCommit("abc")).toBeNull();
  });

  it("commits a trailing-dot draft (Number('1.') === 1) so typing can continue", () => {
    expect(parseCommit("1.")).toBe(1);
  });
});
