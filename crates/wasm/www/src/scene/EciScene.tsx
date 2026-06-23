import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls } from "@react-three/drei";
import type { ChiefGeometry } from "../wasm";
import { scaleAll, type V3 } from "./vec";
import { Arrow } from "./Arrow";

const EARTH_RADIUS_M = 6.378e6;

export function EciScene({ g }: { g: ChiefGeometry }) {
  const k = 1 / g.a; // metres → scene units (a ≈ 1)
  const orbit = scaleAll(g.orbit_eci as V3[], k);
  const arc = g.perigee_arc_eci ? scaleAll(g.perigee_arc_eci as V3[], k) : null;
  const earthR = EARTH_RADIUS_M * k;

  return (
    <div className="canvas3d canvas-eci">
      <Canvas camera={{ position: [2.2, 1.4, 2.2], fov: 45, near: 0.01, far: 100 }}>
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 5, 5]} intensity={0.8} />
        {/* Central body */}
        <mesh>
          <sphereGeometry args={[earthR, 32, 32]} />
          <meshStandardMaterial color="#16324a" wireframe />
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
        {/* Primer arrow at the first burn (static) — amber, DISTINCT from thrust.
            Playback (Task 13) sweeps it along the track. */}
        {g.primer_eci.length > 0 && g.chief_track_eci.length > 0 && (
          <Arrow
            origin={[
              g.chief_track_eci[0][0] * k,
              g.chief_track_eci[0][1] * k,
              g.chief_track_eci[0][2] * k,
            ]}
            dir={g.primer_eci[0] as V3}
            length={0.5}
            color="#ffb454"
          />
        )}
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
