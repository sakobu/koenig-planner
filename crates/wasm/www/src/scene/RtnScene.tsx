import { useMemo } from "react";
import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls, Text } from "@react-three/drei";
import type { ChiefGeometry } from "../lib/wasm";
import { maxRadius, rtnToView, scaleAll, type V3 } from "./vec";
import { RTN_BASIS, RTN_COLORS } from "../lib/rtn";
import { ARROW, SCENE } from "./palette";
import { Arrow } from "./Arrow";

export function RtnScene({ g, sampleIndex }: { g: ChiefGeometry; sampleIndex: number }) {
  // Two curves share one auto-fit: the true transfer (green — the deputy's
  // actual controlled path from its chief-coincident start, kinked at each
  // burn) and the target relative orbit (gray ghost, anchored at t_f, where
  // the transfer lands). A non-zero δa gives open, drifting spirals rather
  // than closed loops: that secular along-track drift is real physics, shown
  // honestly rather than hidden.
  const transfer = g.transfer_track_rtn;
  const target = g.target_track_rtn;
  const rmax = Math.max(1e-6, maxRadius(transfer), maxRadius(target));
  const k = 1 / rmax; // auto-fit meters → ~unit scene
  // Orient with the conventional radial-up / transverse-right / normal-depth
  // axes (see rtnToView), viewed obliquely so the genuinely 3D shape reads
  // honestly. Data stays [radial, transverse, normal]; only the mapping changes.
  const transferCurve = useMemo(() => scaleAll(transfer.map(rtnToView), k), [transfer, k]);
  const targetCurve = useMemo(() => scaleAll(target.map(rtnToView), k), [target, k]);
  const axis = 0.8; // reference-gnomon length; kept short so labels stay inside the viewport

  // Deputy glyph rides the transfer when it is drawable, else the target ghost
  // (degraded transfer), else hides.
  const glyphTrack = transfer.length > 0 ? transfer : target;
  let deputyPos: V3 | null = null;
  if (glyphTrack.length > 0) {
    const v = rtnToView(glyphTrack[sampleIndex] ?? glyphTrack[0]);
    deputyPos = [v[0] * k, v[1] * k, v[2] * k];
  }

  return (
    <div className="canvas3d canvas-rtn">
      <Canvas camera={{ position: [2.0, 1.4, 2.2], fov: 45, near: 0.01, far: 100 }}>
        {/* Lift ambient slightly for the darker console ground. */}
        <ambientLight intensity={0.75} />
        {/* Chief at origin */}
        <mesh>
          <sphereGeometry args={[0.03, 16, 16]} />
          <meshStandardMaterial color={SCENE.spacecraft} />
        </mesh>
        {/* RTN axis gnomon, derived from rtnToView(basis) so the drawn axes and
            labels can never drift from the data's view mapping: T transverse
            (+X), R radial (+Y), N normal (−Z). Physical color binding — R radial
            red, T transverse cyan, N normal amber (see ../rtn). */}
        {(["R", "T", "N"] as const).map((comp) => {
          const dir = rtnToView(RTN_BASIS[comp]);
          const end: V3 = [dir[0] * axis, dir[1] * axis, dir[2] * axis];
          const lpos: V3 = [dir[0] * (axis + 0.1), dir[1] * (axis + 0.1), dir[2] * (axis + 0.1)];
          return (
            <group key={comp}>
              <Line points={[[0, 0, 0], end]} color={RTN_COLORS[comp]} lineWidth={1.5} />
              <Text position={lpos} fontSize={0.12} color={RTN_COLORS[comp]}>
                {comp}
              </Text>
            </group>
          );
        })}
        {/* Target relative orbit — the destination, drawn faint under the transfer. */}
        {targetCurve.length > 1 && (
          <Line points={targetCurve} color={SCENE.targetOrbit} lineWidth={1} />
        )}
        {/* True transfer trajectory — the primary curve; kinks are the burns. */}
        {transferCurve.length > 1 && (
          <Line points={transferCurve} color={SCENE.deputy} lineWidth={2} />
        )}
        {/* Burn nodes + Δv (thrust) arrows — cyan, the same Δv/thrust channel as
            the ECI scene (and the timeline stems). Nodes sit on the transfer's
            kinks (the burn's grid sample); arrows show DIRECTION only (fixed
            length) — per-burn magnitude is read from the Δv-component bars
            (RtnComponents). Both pass through rtnToView so they align with the
            gnomon and the curves; dv_rtn is already the native RTN frame. */}
        {g.maneuver_rtn.map((m, j) => {
          const p = rtnToView(m.position_rtn);
          const pos: V3 = [p[0] * k, p[1] * k, p[2] * k];
          return (
            <group key={j}>
              <mesh position={pos}>
                <sphereGeometry args={[0.03, 12, 12]} />
                <meshStandardMaterial color={SCENE.burn} />
              </mesh>
              <Arrow origin={pos} dir={rtnToView(m.dv_rtn)} length={ARROW.burn} color={SCENE.burn} />
            </group>
          );
        })}
        {/* Deputy glyph + swept primer arrow (amber, like the ECI primer),
            synced to the playback scrubber. At a burn sample the primer arrow
            aligns with that burn's Δv arrow while |p| = 1 — the optimality
            condition made visible. */}
        {deputyPos && (
          <group>
            <mesh position={deputyPos}>
              <sphereGeometry args={[0.04, 16, 16]} />
              <meshStandardMaterial color={SCENE.deputy} />
            </mesh>
            {g.primer_rtn.length > 0 && (
              <Arrow
                origin={deputyPos}
                dir={rtnToView(g.primer_rtn[sampleIndex] ?? g.primer_rtn[0])}
                length={ARROW.primer}
                color={SCENE.primer}
              />
            )}
          </group>
        )}
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
