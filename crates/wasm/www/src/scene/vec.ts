export type V3 = [number, number, number];

export function scaleAll(pts: V3[], k: number): V3[] {
  return pts.map((p) => [p[0] * k, p[1] * k, p[2] * k]);
}

export function maxRadius(pts: V3[]): number {
  let m = 0;
  for (const p of pts) m = Math.max(m, Math.hypot(p[0], p[1], p[2]));
  return m;
}
