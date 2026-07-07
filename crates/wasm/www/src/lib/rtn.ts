/** The R/T/N (radial / transverse / normal) axis triad — the single palette
 *  shared by the component charts and the RTN scene gnomon. Physical binding:
 *  radial red, transverse cyan, normal amber. */
export const RTN_COLORS = { R: "#ff6b6b", T: "#4dd2ff", N: "#ffb454" } as const;
export const RTN_NAME = { R: "radial", T: "transverse", N: "normal" } as const;

/** Unit basis vectors in `[radial, transverse, normal]` order. Map these
 *  through `rtnToView` to place the scene gnomon, so the drawn axes and labels
 *  can never drift from the data's view mapping. */
export const RTN_BASIS: Record<"R" | "T" | "N", [number, number, number]> = {
  R: [1, 0, 0],
  T: [0, 1, 0],
  N: [0, 0, 1],
};
