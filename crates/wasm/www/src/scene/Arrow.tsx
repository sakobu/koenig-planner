import { Line } from "@react-three/drei";
import { Quaternion, Vector3 } from "three";
import type { V3 } from "./vec";

export function Arrow({
  origin,
  dir,
  length,
  color,
}: {
  origin: V3;
  dir: V3;
  length: number;
  color: string;
}) {
  const d = new Vector3(dir[0], dir[1], dir[2]);
  if (d.lengthSq() === 0) return null;
  d.normalize();
  const o = new Vector3(origin[0], origin[1], origin[2]);
  const tip = o.clone().add(d.clone().multiplyScalar(length));
  // Cone default points +Y; rotate it onto d.
  const q = new Quaternion().setFromUnitVectors(new Vector3(0, 1, 0), d);
  const headLen = length * 0.25;
  const coneCenter = tip.clone().sub(d.clone().multiplyScalar(headLen / 2));

  return (
    <group>
      <Line points={[[o.x, o.y, o.z], [tip.x, tip.y, tip.z]]} color={color} lineWidth={2} />
      <mesh position={[coneCenter.x, coneCenter.y, coneCenter.z]} quaternion={q}>
        <coneGeometry args={[length * 0.08, headLen, 12]} />
        <meshStandardMaterial color={color} />
      </mesh>
    </group>
  );
}
