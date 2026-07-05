import { RTN_COLORS, RTN_NAME } from "../rtn";

/** The R/T/N color-legend row shared by the component charts: three swatches
 *  starting at `x`, evenly spread so the last sits 120 units inside the right
 *  edge (`width - padR`), with the row baseline at `y`. */
export function RtnLegend({
  x,
  y,
  width,
  padR,
}: {
  x: number;
  y: number;
  width: number;
  padR: number;
}) {
  const lstep = (width - padR - 120 - x) / 2;
  return (
    <>
      {(["R", "T", "N"] as const).map((comp, k) => {
        const lx = x + k * lstep;
        return (
          <g key={comp}>
            <rect x={lx} y={y - 30} width={11} height={11} rx={2} fill={RTN_COLORS[comp]} />
            <text x={lx + 17} y={y - 20} className="legend-label">
              {RTN_NAME[comp]}
            </text>
          </g>
        );
      })}
    </>
  );
}
