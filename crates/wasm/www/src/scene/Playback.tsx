import { useEffect, useState } from "react";

export function Playback({
  count,
  index,
  setIndex,
}: {
  count: number;
  index: number;
  setIndex: (i: number) => void;
}) {
  const [playing, setPlaying] = useState(false);
  useEffect(() => {
    if (!playing || count <= 1) return;
    const id = setInterval(() => setIndex((index + 1) % count), 60);
    return () => clearInterval(id);
  }, [playing, index, count, setIndex]);

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
    </div>
  );
}
