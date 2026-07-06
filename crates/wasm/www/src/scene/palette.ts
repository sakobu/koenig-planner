/** Semantic body / marker / arrow colors for the 3D scenes. The R/T/N axis
 *  triad lives in `../rtn` (shared with the charts); these are the scene-only
 *  hues, named here so the two scenes can't disagree on a shared one (the
 *  spacecraft marker and primer arrow appear in both). */
export const SCENE = {
  spacecraft: "#dce6f0", // chief / spacecraft marker (both scenes)
  primer: "#ffb454", // swept primer arrow (both scenes) — the amber/primer channel
  burn: "#5cc8ff", // Δv (thrust) nodes + arrows (both scenes) — the cyan Δv/thrust channel
  chiefOrbit: "#7c8b9a", // ECI chief-orbit curve
  perigeeArc: "#ffb454", // ECI perigee attitude-constraint arc
  earthCore: "#0d2336", // ECI central-body core
  earthWire: "#3f86b3", // ECI central-body wireframe
  earthAtmo: "#5cc8ff", // ECI atmosphere rim shell
  deputy: "#5ef2a8", // RTN deputy marker + true-transfer curve
  targetOrbit: "#7c8b9a", // RTN target-orbit ghost (same neutral as the ECI chief orbit)
} as const;

/** Direction-only arrow lengths (scene units; both scenes normalize to a ~unit
 *  view). Kept small so the glyphs read as direction indicators without
 *  dominating the viewport, and shared so the Δv and primer arrows are identical
 *  in size across the ECI and RTN scenes. */
export const ARROW = { burn: 0.22, primer: 0.28 } as const;
