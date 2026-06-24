import type { SolveRequest, CostSpec } from "../wasm";

export type ChiefKey = "a" | "e" | "i" | "raan" | "argp" | "mean_anom";

export function setChief(req: SolveRequest, key: ChiefKey, v: number): SolveRequest {
  return { ...req, chief: { ...req.chief, [key]: v } };
}

export function setWindow(
  req: SolveRequest,
  key: "t_i" | "t_f" | "dt",
  v: number,
): SolveRequest {
  return { ...req, [key]: v };
}

export function setW(req: SolveRequest, idx: number, v: number): SolveRequest {
  const w_metres = [...req.w_metres] as SolveRequest["w_metres"];
  w_metres[idx] = v;
  return { ...req, w_metres };
}

export function setCostType(req: SolveRequest, type: CostSpec["type"]): SolveRequest {
  if (type === "norm2") return { ...req, cost: { type: "norm2" } };
  if (type === "facemax") return { ...req, cost: { type: "facemax" } };
  return { ...req, cost: { type: "piecewise" } };
}

export function setPiecewise(
  req: SolveRequest,
  key: "period" | "t_perigee0",
  v: number | undefined,
): SolveRequest {
  if (req.cost.type !== "piecewise") return req;
  return { ...req, cost: { ...req.cost, [key]: v } };
}
