import { useEffect, useState } from "react";

/** Parse a draft input string to a committable finite number, or null when the
 *  draft is an incomplete entry — empty, a lone "-", a bare "1e" — that should
 *  stay editable without committing, so clearing the field doesn't snap the
 *  value to 0 and a leading minus isn't clobbered to NaN. */
export function parseCommit(draft: string): number | null {
  if (draft.trim() === "") return null;
  const n = Number(draft);
  return Number.isFinite(n) ? n : null;
}

export function NumberField({
  label,
  value,
  onChange,
  min,
  max,
  step,
}: {
  label: string;
  value: number;
  onChange: (v: number) => void;
  min: number;
  max: number;
  step: number;
}) {
  // The number input edits a string draft so intermediate states ("", "-", "1.")
  // stay on screen without committing; only a finite parse reaches onChange.
  const [draft, setDraft] = useState(String(value));
  // Resync when value changes from outside (the slider, a preset), but leave an
  // in-progress edit that already equals value (e.g. "1." while value is 1)
  // alone. Keyed only on value: reacting to draft would fight the user's typing.
  useEffect(() => {
    if (parseCommit(draft) !== value) setDraft(String(value));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [value]);

  const onDraft = (next: string) => {
    setDraft(next);
    const n = parseCommit(next);
    if (n !== null) onChange(n);
  };

  return (
    <div className="field">
      <label>
        {label}
        <input
          type="number"
          step="any"
          value={draft}
          onChange={(e) => onDraft(e.target.value)}
        />
      </label>
      <input
        className="slider"
        type="range"
        aria-label={label}
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
      />
    </div>
  );
}
