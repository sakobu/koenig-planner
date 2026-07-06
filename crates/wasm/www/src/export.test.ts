import { describe, it, expect } from "vitest";
import { toBurnCsv, toPlanJson } from "./export";
import type { SolveRequest, SolveResponse } from "./wasm";

// Minimal structural fixtures — the serializers read only maneuvers + geometry.maneuver_nu.
function respWith(
  maneuvers: { t: number; dv: [number, number, number] }[],
  nu: number[],
): SolveResponse {
  return { maneuvers, geometry: { maneuver_nu: nu } } as unknown as SolveResponse;
}

describe("toBurnCsv", () => {
  it("emits the header plus one full-precision row per maneuver", () => {
    const r = respWith(
      [
        { t: 100, dv: [3, 4, 0] },
        { t: 250, dv: [0.123456789012345, -1, 2] },
      ],
      [0.5, 1.25],
    );
    expect(toBurnCsv(r)).toBe(
      [
        "t_s,dv_R,dv_T,dv_N,dv_mag,nu_rad",
        "100,3,4,0,5,0.5",
        `250,0.123456789012345,-1,2,${Math.hypot(0.123456789012345, -1, 2)},1.25`,
      ].join("\n"),
    );
  });

  it("emits only the header when there are no maneuvers", () => {
    expect(toBurnCsv(respWith([], []))).toBe("t_s,dv_R,dv_T,dv_N,dv_mag,nu_rad");
  });
});

describe("toPlanJson", () => {
  it("round-trips the request and response through JSON", () => {
    const req = { t_i: 0, t_f: 100, dt: 30 } as unknown as SolveRequest;
    const r = respWith([{ t: 10, dv: [1, 2, 3] }], [0.1]);
    expect(JSON.parse(toPlanJson(req, r))).toEqual({ request: req, response: r });
  });
});
