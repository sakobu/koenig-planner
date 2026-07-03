import { describe, it, expect } from "vitest";
import { pickDisplay } from "./outcomeDisplay";
import type { SolveOutcome, SolveResponse } from "./wasm";

const R1 = { total_dv: 1 } as unknown as SolveResponse;
const R2 = { total_dv: 2 } as unknown as SolveResponse;
const OK1 = { status: "ok", value: R1 } as unknown as SolveOutcome;
const ERR = {
  status: "err",
  error: { kind: "bad_request", message: "x" },
} as unknown as SolveOutcome;

describe("pickDisplay", () => {
  it("shows empty when there is no outcome", () => {
    expect(pickDisplay(null, null)).toEqual({ view: "empty" });
  });

  it("shows an ok outcome with no error banner", () => {
    expect(pickDisplay(OK1, null)).toEqual({ view: "ok", r: R1, error: null });
  });

  it("keeps the last-good response and overlays a transient error", () => {
    expect(pickDisplay(ERR, R2)).toEqual({
      view: "ok",
      r: R2,
      error: { kind: "bad_request", message: "x" },
    });
  });

  it("falls back to a bare error when no prior good solve exists", () => {
    expect(pickDisplay(ERR, null)).toEqual({
      view: "error",
      error: { kind: "bad_request", message: "x" },
    });
  });
});
