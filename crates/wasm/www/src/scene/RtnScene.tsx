import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls, Text } from "@react-three/drei";
import type { ChiefGeometry } from "../wasm";
import { maxRadius, scaleAll, type V3 } from "./vec";

export function RtnScene({ g }: { g: ChiefGeometry }) {
  const rel = g.relative_trajectory_rtn as V3[];
  const rmax = Math.max(1e-6, maxRadius(rel));
  const k = 1 / rmax; // auto-fit metres → ~unit scene
  const curve = scaleAll(rel, k);
  const axis = 1.2;

  return (
    <div className="canvas3d canvas-rtn">
      <Canvas camera={{ position: [2, 1.4, 2], fov: 45, near: 0.01, far: 100 }}>
        {/* Lift ambient slightly for the darker console ground. */}
        <ambientLight intensity={0.75} />
        {/* Chief at origin */}
        <mesh>
          <sphereGeometry args={[0.03, 16, 16]} />
          <meshStandardMaterial color="#dce6f0" />
        </mesh>
        {/* RTN axes (R radial, T along-track, N normal) */}
        <Line points={[[0, 0, 0], [axis, 0, 0]]} color="#ff6b6b" lineWidth={1.5} />
        <Line points={[[0, 0, 0], [0, axis, 0]]} color="#4dd2ff" lineWidth={1.5} />
        <Line points={[[0, 0, 0], [0, 0, axis]]} color="#ffb454" lineWidth={1.5} />
        <Text position={[axis + 0.1, 0, 0]} fontSize={0.12} color="#ff6b6b">R</Text>
        <Text position={[0, axis + 0.1, 0]} fontSize={0.12} color="#4dd2ff">T</Text>
        <Text position={[0, 0, axis + 0.1]} fontSize={0.12} color="#ffb454">N</Text>
        {/* Deputy relative orbit */}
        <Line points={curve} color="#5ef2a8" lineWidth={2} />
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
