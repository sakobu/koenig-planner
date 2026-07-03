import { useEffect, useState, type Dispatch, type SetStateAction } from "react";

// Fixed visual tick; the speed multiplier advances that many samples per tick
// (rather than shrinking the interval) so high speeds stay smooth instead of
// hitting the browser's minimum-timer clamp.
const TICK_MS = 50;
const SPEEDS = [1, 8, 64] as const;

export function Playback({
  count,
  index,
  setIndex,
}: {
  count: number;
  index: number;
  setIndex: Dispatch<SetStateAction<number>>;
}) {
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(1);
  useEffect(() => {
    if (!playing || count <= 1) return;
    const id = setInterval(() => setIndex((i) => (i + speed) % count), TICK_MS);
    return () => clearInterval(id);
  }, [playing, count, setIndex, speed]);

  return (
    <div className="playback">
      <button type="button" onClick={() => setPlaying((p) => !p)}>
        {playing ? "❚❚" : "▶"}
      </button>
      <input
        type="range"
        min={0}
        max={Math.max(0, count - 1)}
        step={1}
        value={index}
        onChange={(e) => {
          setPlaying(false);
          setIndex(Number(e.target.value));
        }}
      />
      <select
        className="speed"
        aria-label="playback speed"
        value={speed}
        onChange={(e) => setSpeed(Number(e.target.value))}
      >
        {SPEEDS.map((s) => (
          <option key={s} value={s}>{`${s}×`}</option>
        ))}
      </select>
    </div>
  );
}
