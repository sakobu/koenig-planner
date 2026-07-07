export type V3 = [number, number, number];

export function scaleAll(pts: V3[], k: number): V3[] {
  return pts.map((p) => [p[0] * k, p[1] * k, p[2] * k]);
}

export function maxRadius(pts: V3[]): number {
  let m = 0;
  for (const p of pts) m = Math.max(m, Math.hypot(p[0], p[1], p[2]));
  return m;
}

/** Map an RTN component triple `[radial, transverse, normal]` into the scene's
 *  view frame: transverse → X (horizontal), radial → Y (vertical), normal → −Z
 *  (depth). A proper rotation (det +1), so the relative orbit is rotated, not
 *  mirrored — the deputy's sense of motion is preserved. With radial up and the
 *  transverse (along-track) axis to the right, an in-plane-dominated orbit reads
 *  as the conventional tilted 2:1 ellipse; a cross-track-dominated one reads
 *  honestly as a 3D loop. */
export function rtnToView(v: V3): V3 {
  return [v[1], v[0], -v[2]];
}

/** Map an ECI vector `[x, y, z]` (Earth-centred inertial, Z = north pole, i.e.
 *  Z-up) into the scene's three.js view frame (Y-up): ECI Z → view Y (up). A
 *  proper rotation (det +1, a −90° turn about X), so orbit sense/handedness is
 *  preserved. With the pole vertical, an equatorial orbit reads as a horizontal
 *  ring and a polar one as a vertical loop — matching how these orbits are
 *  conventionally pictured. (three.js is Y-up; feeding ECI in raw would lay the
 *  pole sideways and stand the equatorial plane upright.) */
export function eciToView(v: V3): V3 {
  return [v[0], v[2], -v[1]];
}
