import { useMemo } from "react";
import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls, Text } from "@react-three/drei";
import type { ChiefGeometry } from "../wasm";
import { maxRadius, rtnToView, scaleAll, type V3 } from "./vec";
import { Arrow } from "./Arrow";

export function RtnScene({ g, sampleIndex }: { g: ChiefGeometry; sampleIndex: number }) {
  // The deputy track is sampled on the playback grid over the FULL mission
  // window (several chief periods), and is BOTH the drawn curve and the glyph
  // source — so the glyph rides the line exactly for the entire scrub. A
  // non-zero δa gives the deputy a slightly different period, so the curve is an
  // open, drifting spiral rather than a single closed loop: that secular
  // along-track drift is real physics, shown honestly rather than hidden.
  const track = g.deputy_track_rtn as V3[];
  const rmax = Math.max(1e-6, maxRadius(track)); // rotation-invariant, so map order is irrelevant
  const k = 1 / rmax; // auto-fit meters → ~unit scene
  // Orient with the conventional radial-up / transverse-right / normal-depth
  // axes (see rtnToView), viewed obliquely so the genuinely 3D shape reads
  // honestly: an in-plane-dominated orbit shows the tilted 2:1 ellipse, a
  // cross-track-dominated one (e.g. the paper's δi-heavy example) reads as a 3D
  // loop. Data stays [radial, transverse, normal]; only the mapping changes.
  const curve = useMemo(() => scaleAll(track.map(rtnToView), k), [track, k]);
  const axis = 0.8; // reference-gnomon length; kept short so labels stay inside the viewport

  // Deputy glyph: position at the current playback sample, same scale/mapping as the curve.
  const clampedIdx = Math.min(sampleIndex, Math.max(0, track.length - 1));
  let deputyPos: V3 | null = null;
  if (track.length > 0) {
    const v = rtnToView(track[clampedIdx]);
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
          <meshStandardMaterial color="#dce6f0" />
        </mesh>
        {/* RTN axes, reoriented: T transverse (horizontal, +X), R radial
            (vertical, +Y), N normal (depth, −Z). Colors keep their physical
            binding: R radial red, T transverse cyan, N normal amber. */}
        <Line points={[[0, 0, 0], [axis, 0, 0]]} color="#4dd2ff" lineWidth={1.5} />
        <Line points={[[0, 0, 0], [0, axis, 0]]} color="#ff6b6b" lineWidth={1.5} />
        <Line points={[[0, 0, 0], [0, 0, -axis]]} color="#ffb454" lineWidth={1.5} />
        <Text position={[axis + 0.1, 0, 0]} fontSize={0.12} color="#4dd2ff">T</Text>
        <Text position={[0, axis + 0.1, 0]} fontSize={0.12} color="#ff6b6b">R</Text>
        <Text position={[0, 0, -axis - 0.1]} fontSize={0.12} color="#ffb454">N</Text>
        {/* Deputy relative orbit */}
        <Line points={curve} color="#5ef2a8" lineWidth={2} />
        {/* Burn nodes + Δv (thrust) arrows — violet. Arrows show DIRECTION only
            (fixed length); per-burn magnitude is read from the Δv-component bars
            (RtnComponents), matching the ECI scene. Both the position and the Δv
            pass through rtnToView so they align with the gnomon and the deputy
            curve; dv_rtn is already the native RTN frame, so no extra rotation.
            The node sits on the deputy curve as a schematic anchor (see
            geometry.rs) — only the arrow direction is exact. */}
        {g.maneuver_rtn.map((m, j) => {
          const p = rtnToView(m.position_rtn as V3);
          const pos: V3 = [p[0] * k, p[1] * k, p[2] * k];
          return (
            <group key={j}>
              <mesh position={pos}>
                <sphereGeometry args={[0.03, 12, 12]} />
                <meshStandardMaterial color="#c792ff" />
              </mesh>
              <Arrow origin={pos} dir={rtnToView(m.dv_rtn as V3)} length={0.3} color="#c792ff" />
            </group>
          );
        })}
        {/* Deputy glyph + swept primer arrow (amber, like the ECI primer),
            synced to the playback scrubber. */}
        {deputyPos && (
          <group>
            <mesh position={deputyPos}>
              <sphereGeometry args={[0.04, 16, 16]} />
              <meshStandardMaterial color="#5ef2a8" />
            </mesh>
            {g.primer_rtn.length > 0 && (
              <Arrow
                origin={deputyPos}
                dir={rtnToView((g.primer_rtn[clampedIdx] ?? g.primer_rtn[0]) as V3)}
                length={0.4}
                color="#ffb454"
              />
            )}
          </group>
        )}
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
