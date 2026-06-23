import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls, Stars } from "@react-three/drei";
import type { ChiefGeometry } from "../wasm";
import { scaleAll, type V3 } from "./vec";
import { Arrow } from "./Arrow";

const EARTH_RADIUS_M = 6.378e6;

export function EciScene({ g, sampleIndex }: { g: ChiefGeometry; sampleIndex: number }) {
  const k = 1 / g.a; // metres → scene units (a ≈ 1)
  const orbit = scaleAll(g.orbit_eci as V3[], k);
  const arc = g.perigee_arc_eci ? scaleAll(g.perigee_arc_eci as V3[], k) : null;
  const earthR = EARTH_RADIUS_M * k;

  return (
    <div className="canvas3d canvas-eci">
      <Canvas camera={{ position: [2.2, 1.4, 2.2], fov: 45, near: 0.01, far: 100 }}>
        <ambientLight intensity={0.55} />
        <directionalLight position={[5, 5, 5]} intensity={0.8} />
        {/* Faint static starfield — instrument backdrop, not noise. */}
        <Stars radius={6} depth={8} count={1200} factor={0.15} saturation={0} fade speed={0} />
        {/* Central body — schematic deep-instrument-blue wireframe. */}
        <mesh>
          <sphereGeometry args={[earthR, 32, 32]} />
          <meshStandardMaterial color="#123a52" wireframe />
        </mesh>
        {/* ECI reference axes */}
        <axesHelper args={[1.6]} />
        {/* Chief orbit */}
        <Line points={orbit} color="#7c8b9a" lineWidth={1.5} />
        {/* FaceMax perigee-window arc (piecewise only) */}
        {arc && <Line points={arc} color="#ffb454" lineWidth={3} />}
        {/* Burn nodes + Δv (thrust) arrows — cyan. */}
        {g.maneuver_eci.map((m, j) => {
          const pos: V3 = [m.position_eci[0] * k, m.position_eci[1] * k, m.position_eci[2] * k];
          return (
            <group key={j}>
              <mesh position={pos}>
                <sphereGeometry args={[0.02, 12, 12]} />
                <meshStandardMaterial color="#5cc8ff" />
              </mesh>
              <Arrow origin={pos} dir={m.dv_eci as V3} length={0.35} color="#5cc8ff" />
            </group>
          );
        })}
        {/* Spacecraft + swept primer at the current playback sample. */}
        {g.chief_track_eci.length > 0 &&
          (() => {
            const i = Math.min(sampleIndex, g.chief_track_eci.length - 1);
            const c = g.chief_track_eci[i];
            const pos: V3 = [c[0] * k, c[1] * k, c[2] * k];
            const primer = g.primer_eci[i] ?? g.primer_eci[0];
            return (
              <group>
                <mesh position={pos}>
                  <sphereGeometry args={[0.03, 16, 16]} />
                  <meshStandardMaterial color="#dce6f0" />
                </mesh>
                {primer && <Arrow origin={pos} dir={primer as V3} length={0.5} color="#ffb454" />}
              </group>
            );
          })()}
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
