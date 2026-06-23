import { Canvas } from "@react-three/fiber";
import { Line, OrbitControls, Text } from "@react-three/drei";
import type { ChiefGeometry } from "../wasm";
import { maxRadius, scaleAll, type V3 } from "./vec";

export function RtnScene({ g, sampleIndex }: { g: ChiefGeometry; sampleIndex: number }) {
  const rel = g.relative_trajectory_rtn as V3[];
  const rmax = Math.max(1e-6, maxRadius(rel));
  const k = 1 / rmax; // auto-fit metres → ~unit scene
  const curve = scaleAll(rel, k);
  const axis = 1.2;

  // Deputy glyph: position at the current playback sample, same scale as the orbit loop.
  const track = g.deputy_track_rtn as V3[];
  const clampedIdx = Math.min(sampleIndex, Math.max(0, track.length - 1));
  const deputyPos: V3 | null =
    track.length > 0
      ? [track[clampedIdx][0] * k, track[clampedIdx][1] * k, track[clampedIdx][2] * k]
      : null;

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
        {/* Deputy glyph synced to the playback scrubber */}
        {deputyPos && (
          <mesh position={deputyPos}>
            <sphereGeometry args={[0.04, 16, 16]} />
            <meshStandardMaterial color="#5ef2a8" />
          </mesh>
        )}
        <OrbitControls enablePan enableZoom enableRotate />
      </Canvas>
    </div>
  );
}
