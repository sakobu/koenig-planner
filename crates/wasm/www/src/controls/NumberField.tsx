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

/** Commit policy for an OPTIONAL numeric field: an empty draft commits
 *  `undefined` ("auto"); a finite parse commits that number; an incomplete
 *  draft ("-", "1e") does not commit, so it stays editable and never reaches the
 *  solver as NaN. */
export function optionalCommit(
  draft: string,
): { commit: false } | { commit: true; value: number | undefined } {
  if (draft.trim() === "") return { commit: true, value: undefined };
  const n = parseCommit(draft);
  return n === null ? { commit: false } : { commit: true, value: n };
}

/** Owns the string-draft mirror of a numeric field so intermediate states ("",
 *  "-", "1.") stay on screen and only a finite parse commits. Resyncs when the
 *  value changes from outside (the slider, a preset) but leaves an in-progress
 *  edit that already equals the value alone. Keyed only on `value`: reacting to
 *  the draft would fight the user's typing. `undefined` shows as an empty draft
 *  (the optional/"auto" case). */
export function useNumberDraft(
  value: number | undefined,
): readonly [string, (next: string) => void] {
  const [draft, setDraft] = useState(value === undefined ? "" : String(value));
  useEffect(() => {
    if (parseCommit(draft) !== (value ?? null)) {
      setDraft(value === undefined ? "" : String(value));
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [value]);
  return [draft, setDraft] as const;
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
  const [draft, setDraft] = useNumberDraft(value);

  const onDraft = (next: string) => {
    setDraft(next);
    const n = parseCommit(next);
    if (n !== null) onChange(n); // live per-keystroke so the 3D tracks typing
  };
  // Clamp to the slider's [min,max] on blur (not per-keystroke, which would make
  // a min-bounded field impossible to type into), so the two editors of one value
  // agree on range instead of the text input silently admitting a forbidden value.
  const onBlur = () => {
    const n = parseCommit(draft);
    if (n === null) return;
    const clamped = Math.min(max, Math.max(min, n));
    if (clamped !== n) {
      setDraft(String(clamped));
      onChange(clamped);
    }
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
          onBlur={onBlur}
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

/** A nullable scalar field (no slider, no range): empty commits `undefined`
 *  ("auto"), a finite parse commits the number. Shares `useNumberDraft` so it
 *  gets the same NaN-proof draft buffering as {@link NumberField} — the piecewise
 *  cost's optional `period` / `t_perigee0` inputs route through here rather than
 *  a second, un-hardened `Number()` parse. */
export function OptionalNumberField({
  label,
  value,
  onChange,
  placeholder = "auto",
}: {
  label: string;
  value: number | undefined;
  onChange: (v: number | undefined) => void;
  placeholder?: string;
}) {
  const [draft, setDraft] = useNumberDraft(value);

  const onDraft = (next: string) => {
    setDraft(next);
    const c = optionalCommit(next);
    if (c.commit) onChange(c.value);
  };

  return (
    <label>
      {label}
      <input
        type="number"
        step="any"
        placeholder={placeholder}
        value={draft}
        onChange={(e) => onDraft(e.target.value)}
      />
    </label>
  );
}
