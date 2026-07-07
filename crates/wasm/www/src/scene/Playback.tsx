import { useEffect, useState, type CSSProperties, type Dispatch, type SetStateAction } from "react";
import {
  burnTickFractions,
  fmtHours,
  fmtNu,
  fmtOrbit,
  nextBurnFrame,
  prevBurnFrame,
} from "./playbackUtil";

// Fixed visual tick; the speed multiplier advances that many samples per tick
// (rather than shrinking the interval) so high speeds stay smooth instead of
// hitting the browser's minimum-timer clamp.
const TICK_MS = 50;
const SPEEDS = [1, 8, 64] as const;

export function Playback({
  count,
  index,
  setIndex,
  times,
  nu,
  burnTimes,
  period,
}: {
  count: number;
  index: number;
  setIndex: Dispatch<SetStateAction<number>>;
  /** Absolute grid times (primer_times); times[k] pairs with sample k. */
  times: number[];
  /** Chief true anomaly per sample [rad] (chief_nu_track). */
  nu: number[];
  /** Absolute burn times (maneuvers[].t). */
  burnTimes: number[];
  /** Chief period [s] for the orbit counter. */
  period: number;
}) {
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(1);
  useEffect(() => {
    if (!playing || count <= 1) return;
    const id = setInterval(() => setIndex((i) => (i + speed) % count), TICK_MS);
    return () => clearInterval(id);
  }, [playing, count, setIndex, speed]);

  const frame = Math.min(index, Math.max(0, count - 1));
  const maxIdx = Math.max(1, count - 1);
  const fillFrac = frame / maxIdx;
  const fillPct = (fillFrac * 100).toFixed(3);
  // The drawn thumb and the burn ticks share this exact center-travel expression,
  // so a step-to-burn lands the thumb on the tick (see style.css --pb-thumb).
  const atFrac = (f: number) => `calc(${f.toFixed(5)} * (100% - var(--pb-thumb)) + var(--pb-thumb) / 2)`;
  const t = times.length ? times[frame] ?? times[0] : 0;
  const t0 = times.length ? times[0] : 0;
  const ticks = burnTickFractions(times, burnTimes);
  const prev = prevBurnFrame(times, burnTimes, frame);
  const next = nextBurnFrame(times, burnTimes, frame);
  // Jumping to a burn pauses, like manual scrubbing: the point is to inspect
  // the primer-alignment moment, not fly past it.
  const jump = (to: number | null) => {
    if (to === null) return;
    setPlaying(false);
    setIndex(to);
  };

  return (
    <div className="playback">
      <div className="playback-controls">
        <button type="button" aria-label="previous burn" disabled={prev === null} onClick={() => jump(prev)}>
          ⏮
        </button>
        <button
          type="button"
          aria-label={playing ? "pause" : "play"}
          onClick={() => setPlaying((p) => !p)}
        >
          {playing ? "❚❚" : "▶"}
        </button>
        <button type="button" aria-label="next burn" disabled={next === null} onClick={() => jump(next)}>
          ⏭
        </button>
        <div className="playback-track">
          <input
            type="range"
            aria-label="playback position"
            min={0}
            max={Math.max(0, count - 1)}
            step={1}
            value={index}
            style={{ "--pb-val": `${fillPct}%` } as CSSProperties}
            onChange={(e) => {
              setPlaying(false);
              setIndex(Number(e.target.value));
            }}
          />
          {/* Decorative burn markers on the scrub track: the readout and the
              plan table carry the values, so these stay out of the a11y tree. */}
          <div className="playback-ticks" aria-hidden="true">
            {ticks.map((f, j) => (
              <span key={j} className="playback-tick" style={{ left: atFrac(f) }} />
            ))}
          </div>
          <span className="pb-thumb" aria-hidden="true" style={{ left: atFrac(fillFrac) }} />
        </div>
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
      {times.length > 0 && (
        <div className="playback-readout">
          {fmtHours(t, t0)} · {fmtOrbit(t, t0, period)} · {nu.length > 0 ? fmtNu(nu[frame] ?? nu[0]) : "ν —"}
        </div>
      )}
    </div>
  );
}
