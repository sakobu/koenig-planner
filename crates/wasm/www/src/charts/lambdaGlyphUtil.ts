/** Pure geometry for the λ dual-certificate glyph: one zero-centered signed bar
 *  per ROE component, length ∝ |λ_i| normalized by the largest |λ|. */
import { maxAbs } from "../lib/svgUtil";

/** ROE component labels, in the `w_meters` / `target_roe` / `lambda` order. */
export const ROE_LABELS = ["δa", "δλ", "δe_x", "δe_y", "δi_x", "δi_y"] as const;

export interface LambdaBar {
  label: string;
  /** Raw dual component (shown at the bar tip). */
  value: number;
  /** Signed fraction of the half-width in [−1, 1] (value / max|λ|). */
  frac: number;
}

/** Bars for the 6-vector `lambda`. `frac` is normalized by the largest
 *  |component| (floored, so an all-zero λ yields all-zero fracs, never NaN). */
export function lambdaBars(lambda: number[]): LambdaBar[] {
  const norm = maxAbs(lambda, 1e-12);
  return lambda.map((value, i) => ({
    label: ROE_LABELS[i] ?? `λ${i}`,
    value,
    frac: value / norm,
  }));
}
