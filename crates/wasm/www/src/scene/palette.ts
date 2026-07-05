/** Semantic body / marker / arrow colors for the 3D scenes. The R/T/N axis
 *  triad lives in `../rtn` (shared with the charts); these are the scene-only
 *  hues, named here so the two scenes can't disagree on a shared one (the
 *  spacecraft marker and primer arrow appear in both). */
export const SCENE = {
  spacecraft: "#dce6f0", // chief / spacecraft marker (both scenes)
  primer: "#ffb454", // swept primer arrow (both scenes)
  eciBurn: "#5cc8ff", // ECI Δv nodes + arrows
  chiefOrbit: "#7c8b9a", // ECI chief-orbit curve
  perigeeArc: "#ffb454", // ECI perigee attitude-constraint arc
  earthCore: "#0d2336", // ECI central-body core
  earthWire: "#3f86b3", // ECI central-body wireframe
  earthAtmo: "#5cc8ff", // ECI atmosphere rim shell
  rtnBurn: "#c792ff", // RTN Δv nodes + arrows
  deputy: "#5ef2a8", // RTN deputy marker + relative-orbit curve
} as const;
