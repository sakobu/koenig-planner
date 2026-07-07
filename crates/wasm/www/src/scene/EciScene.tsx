import { useMemo } from "react";
import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls, Stars } from "@react-three/drei";
import { BackSide } from "three";
import type { ChiefGeometry } from "../lib/wasm";
import { eciToView, scaleAll, type V3 } from "./vec";
import { ARROW, SCENE } from "./palette";
import { Arrow } from "./Arrow";

const EARTH_RADIUS_M = 6.378e6;

export function EciScene({ g, sampleIndex }: { g: ChiefGeometry; sampleIndex: number }) {
  const k = 1 / g.a; // meters → scene units (a ≈ 1)
  // Lift every ECI-space quantity into the view frame (pole up) via eciToView.
  const orbit = useMemo(() => scaleAll(g.orbit_eci, k).map(eciToView), [g, k]);
  const arc = useMemo(
    () => (g.perigee_arc_eci ? scaleAll(g.perigee_arc_eci, k).map(eciToView) : null),
    [g, k],
  );
  const earthR = EARTH_RADIUS_M * k;

  return (
    <div className="canvas3d canvas-eci">
      <Canvas camera={{ position: [2.2, 1.4, 2.2], fov: 45, near: 0.01, far: 100 }}>
        <ambientLight intensity={0.55} />
        <directionalLight position={[5, 5, 5]} intensity={0.8} />
        {/* Faint static starfield — instrument backdrop, not noise. */}
        <Stars radius={6} depth={8} count={1200} factor={0.15} saturation={0} fade speed={0} />
        {/* Central body — schematic instrument globe: a solid dark core for
            mass, a lighting-INDEPENDENT steel-cyan wireframe (meshBasic, so the
            far side never goes dark against the near-black ground), and a faint
            cyan atmosphere rim (back-face shell). */}
        <group>
          <mesh>
            <sphereGeometry args={[earthR * 0.99, 32, 32]} />
            <meshStandardMaterial color={SCENE.earthCore} />
          </mesh>
          <mesh>
            <sphereGeometry args={[earthR, 24, 24]} />
            <meshBasicMaterial color={SCENE.earthWire} wireframe />
          </mesh>
          <mesh>
            <sphereGeometry args={[earthR * 1.06, 24, 24]} />
            <meshBasicMaterial color={SCENE.earthAtmo} transparent opacity={0.05} side={BackSide} />
          </mesh>
        </group>
        {/* Reference triad — after the eciToView lift, green (up) is ECI +Z (north pole). */}
        <axesHelper args={[1.6]} />
        {/* Chief orbit */}
        <Line points={orbit} color={SCENE.chiefOrbit} lineWidth={1.5} />
        {/* Perigee attitude-constraint window arc (piecewise only — eq. 49's T1,
            where the cost switches to FaceMax; Norm2 elsewhere) */}
        {arc && <Line points={arc} color={SCENE.perigeeArc} lineWidth={3} />}
        {/* Burn nodes + Δv (thrust) arrows — cyan. Arrows show DIRECTION only
            (fixed length); per-burn magnitude is read from the Δv-component
            bars (RtnComponents). Same for the amber primer arrow below. */}
        {g.maneuver_eci.map((m, j) => {
          const pos: V3 = eciToView([m.position_eci[0] * k, m.position_eci[1] * k, m.position_eci[2] * k]);
          return (
            <group key={j}>
              <mesh position={pos}>
                <sphereGeometry args={[0.02, 12, 12]} />
                <meshStandardMaterial color={SCENE.burn} />
              </mesh>
              <Arrow origin={pos} dir={eciToView(m.dv_eci)} length={ARROW.burn} color={SCENE.burn} />
            </group>
          );
        })}
        {/* Spacecraft + swept primer at the current playback sample. */}
        {g.chief_track_eci.length > 0 &&
          (() => {
            const c = g.chief_track_eci[sampleIndex];
            const pos: V3 = eciToView([c[0] * k, c[1] * k, c[2] * k]);
            const primer = g.primer_eci[sampleIndex] ?? g.primer_eci[0];
            return (
              <group>
                <mesh position={pos}>
                  <sphereGeometry args={[0.03, 16, 16]} />
                  <meshStandardMaterial color={SCENE.spacecraft} />
                </mesh>
                {primer && <Arrow origin={pos} dir={eciToView(primer)} length={ARROW.primer} color={SCENE.primer} />}
              </group>
            );
          })()}
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
