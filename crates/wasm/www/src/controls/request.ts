import type { SolveRequest, CostSpec } from "../lib/wasm";

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
  const w_meters = [...req.w_meters] as SolveRequest["w_meters"];
  w_meters[idx] = v;
  return { ...req, w_meters };
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
