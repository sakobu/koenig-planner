import type {
  SolveOutcome,
  SolveResponse,
  ChiefGeometry,
} from "koenig-damico-planner-wasm";

const NS = "http://www.w3.org/2000/svg";
const RTN_COLORS = { R: "#ff6b6b", T: "#4dd2ff", N: "#ffb454" };

function el(
  tag: string,
  attrs: Record<string, string | number>,
  text?: string,
): SVGElement {
  const node = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, String(v));
  if (text !== undefined) node.textContent = text;
  return node as SVGElement;
}

function svg(w: number, h: number, cls: string): SVGSVGElement {
  return el("svg", {
    viewBox: `0 0 ${w} ${h}`,
    width: "100%",
    preserveAspectRatio: "xMidYMid meet",
    class: cls,
  }) as SVGSVGElement;
}

// Round a raw step up to a 1/2/5 ×10ⁿ "nice" increment for axis ticks.
function niceStep(raw: number): number {
  const exp = Math.floor(Math.log10(raw));
  const f = raw / 10 ** exp;
  const nf = f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10;
  return nf * 10 ** exp;
}

function kpis(r: SolveResponse): HTMLElement {
  const box = document.createElement("div");
  box.className = "kpis";
  const cells: [string, string][] = [
    ["Δv cost", `${r.total_dv.toFixed(4)} m/s`],
    ["maneuvers", String(r.maneuvers.length)],
    ["iterations", String(r.iterations)],
    ["residual", r.residual.toExponential(2)],
  ];
  for (const [label, value] of cells) {
    const c = document.createElement("div");
    c.className = "kpi";
    c.innerHTML = `<span class="k-label">${label}</span><span class="k-value">${value}</span>`;
    box.appendChild(c);
  }
  return box;
}

function timeline(r: SolveResponse, t_i: number, t_f: number): SVGSVGElement {
  // Wide, short aspect — natural for a time series and keeps the panel from
  // ballooning when the column stretches to full width.
  const W = 760,
    H = 300;
  const padL = 58,
    padR = 30,
    padT = 46, // headroom for the stacked "mnvr N" + value labels over a stem
    padB = 44;
  const yBase = H - padB,
    yTop = padT;
  const plotH = yBase - yTop,
    plotW = W - padL - padR;
  const s = svg(W, H, "chart chart-timeline");
  const mags = r.maneuvers.map((m) => Math.hypot(m.dv[0], m.dv[1], m.dv[2]));
  const maxMag = Math.max(1e-12, ...mags);
  // Round the y-axis to human-friendly increments (0.01, 0.02, …) instead of
  // raw fractions of the data max.
  const step = niceStep(maxMag / 4);
  const domainMax = Math.max(step, Math.ceil(maxMag / step) * step);
  const inset = 0.1 * plotW; // keep first/last burns (and their tags) inset
  const span = Math.max(1e-9, t_f - t_i);
  const x = (t: number) =>
    padL + inset + ((t - t_i) / span) * (plotW - 2 * inset);
  const y = (mag: number) => yBase - (mag / domainMax) * plotH;

  // Horizontal gridlines + magnitude ticks (rounded, fixed-decimal).
  for (let v = 0; v <= domainMax + step / 2; v += step) {
    const gy = y(v);
    s.appendChild(
      el("line", {
        x1: padL,
        y1: gy,
        x2: W - padR,
        y2: gy,
        class: v === 0 ? "axis" : "grid",
      }),
    );
    s.appendChild(
      el(
        "text",
        { x: padL - 10, y: gy + 3.5, class: "axis-label", "text-anchor": "end" },
        v.toFixed(4),
      ),
    );
  }
  // Horizontal y-axis caption in the top-left margin — readable, never
  // perpendicular, and clear of the tick column below it.
  s.appendChild(
    el(
      "text",
      { x: 6, y: 15, class: "axis-title", "text-anchor": "start" },
      "|Δv|  [m/s]",
    ),
  );
  // X-axis endpoints (the first/last burn times) + centered title.
  s.appendChild(
    el(
      "text",
      { x: x(t_i), y: yBase + 18, class: "axis-label", "text-anchor": "middle" },
      t_i.toFixed(0),
    ),
  );
  s.appendChild(
    el(
      "text",
      { x: x(t_f), y: yBase + 18, class: "axis-label", "text-anchor": "middle" },
      t_f.toFixed(0),
    ),
  );
  s.appendChild(
    el(
      "text",
      {
        x: padL + plotW / 2,
        y: yBase + 35,
        class: "axis-title",
        "text-anchor": "middle",
      },
      "burn time  [s]",
    ),
  );
  r.maneuvers.forEach((m, j) => {
    const mag = mags[j];
    const mx = x(m.t),
      my = y(mag);
    s.appendChild(
      el("line", { x1: mx, y1: yBase, x2: mx, y2: my, class: "stem" }),
    );
    s.appendChild(el("circle", { cx: mx, cy: my, r: 4, class: "stem-dot" }));
    // Value over each dot, with the maneuver index above it — the same index
    // the R/T/N panel uses, so the two charts read together.
    s.appendChild(
      el(
        "text",
        { x: mx, y: my - 11, class: "stem-label", "text-anchor": "middle" },
        mag.toFixed(4),
      ),
    );
    s.appendChild(
      el(
        "text",
        { x: mx, y: my - 25, class: "mnvr-tag", "text-anchor": "middle" },
        `mnvr ${j + 1}`,
      ),
    );
  });
  return s;
}

const RTN_NAME = { R: "radial", T: "tangential", N: "normal" } as const;

function rtnBars(r: SolveResponse): SVGSVGElement {
  // Diverging bars about a centered zero axis. A left gutter holds the maneuver
  // labels and the value (channel + signed magnitude) rides at each bar tip, so
  // every bar is legible without leaning on the color key alone.
  const n = r.maneuvers.length;
  const W = 760,
    rowH = 54,
    padL = 84,
    padR = 104,
    padT = 42,
    padB = 26;
  const H = padT + n * rowH + padB;
  const s = svg(W, H, "chart chart-rtn");
  const maxComp = Math.max(
    1e-12,
    ...r.maneuvers.flatMap((m) => m.dv.map(Math.abs)),
  );
  const plotW = W - padL - padR;
  const cx = padL + plotW / 2; // shared zero axis
  const labelRoom = 72; // reserve space at the bar tips for the value text
  // Round the half-range so the gridlines land on human-friendly values.
  const tickStep = niceStep(maxComp / 3);
  const domainMax = Math.max(tickStep, Math.ceil(maxComp / tickStep) * tickStep);
  const scale = (plotW / 2 - labelRoom) / domainMax;
  const bh = 11,
    gap = 5,
    blockH = 3 * bh + 2 * gap;
  const axisTop = padT - 6,
    axisBot = padT + n * rowH + 2;

  // Color key (radial / tangential / normal), spread evenly across the width
  // with the unit tucked into the right corner.
  const lstep = (W - padR - 120 - padL) / 2;
  (["R", "T", "N"] as const).forEach((comp, k) => {
    const lx = padL + k * lstep;
    s.appendChild(
      el("rect", {
        x: lx,
        y: padT - 30,
        width: 11,
        height: 11,
        rx: 2,
        fill: RTN_COLORS[comp],
      }),
    );
    s.appendChild(
      el("text", { x: lx + 17, y: padT - 20, class: "legend-label" }, RTN_NAME[comp]),
    );
  });
  s.appendChild(
    el(
      "text",
      { x: W - padR, y: padT - 20, class: "axis-label", "text-anchor": "end" },
      "[m/s]",
    ),
  );

  // Quantitative scale — faint signed gridlines so bar lengths read true.
  for (let v = -domainMax; v <= domainMax + tickStep / 2; v += tickStep) {
    const gx2 = cx + v * scale;
    s.appendChild(
      el("line", {
        x1: gx2,
        y1: axisTop,
        x2: gx2,
        y2: axisBot,
        class: Math.abs(v) < tickStep / 2 ? "zero-axis" : "grid",
      }),
    );
    s.appendChild(
      el(
        "text",
        { x: gx2, y: axisBot + 16, class: "axis-label", "text-anchor": "middle" },
        Math.abs(v) < tickStep / 2 ? "0" : v.toFixed(3),
      ),
    );
  }

  r.maneuvers.forEach((m, j) => {
    const yc = padT + j * rowH + rowH / 2;
    const top = yc - blockH / 2;
    s.appendChild(
      el(
        "text",
        { x: padL - 16, y: yc + 4, class: "row-label", "text-anchor": "end" },
        `mnvr ${j + 1}`,
      ),
    );
    (["R", "T", "N"] as const).forEach((comp, k) => {
      const v = m.dv[k];
      const by = top + k * (bh + gap);
      const len = Math.abs(v) * scale;
      const pos = v >= 0;
      s.appendChild(
        el("rect", {
          x: pos ? cx : cx - len,
          y: by,
          width: Math.max(len, 0.75),
          height: bh,
          rx: 1.5,
          fill: RTN_COLORS[comp],
        }),
      );
      s.appendChild(
        el(
          "text",
          {
            x: pos ? cx + len + 7 : cx - len - 7,
            y: by + bh - 1.5,
            class: "val-label",
            "text-anchor": pos ? "start" : "end",
          },
          `${comp} ${pos ? "+" : "−"}${Math.abs(v).toFixed(4)}`,
        ),
      );
    });
  });
  return s;
}

function orbit(g: ChiefGeometry): SVGSVGElement {
  const W = 320,
    H = 300;
  const s = svg(W, H, "chart chart-orbit");
  const e = g.e,
    a = 1;
  const rOf = (nu: number) => (a * (1 - e * e)) / (1 + e * Math.cos(nu));
  const gx = (nu: number) => rOf(nu) * Math.cos(nu); // focus at the origin
  const gy = (nu: number) => rOf(nu) * Math.sin(nu);

  // An e=0.7 ellipse is wide, short, and offset from its focus, so a
  // focus-centered fit wastes most of the panel. Measure the ellipse's true
  // bounding box and fit *that* — then the figure fills the frame and centers.
  let minX = Infinity,
    maxX = -Infinity,
    minY = Infinity,
    maxY = -Infinity;
  const N = 240;
  for (let i = 0; i <= N; i++) {
    const nu = -Math.PI + (i / N) * 2 * Math.PI;
    const X = gx(nu),
      Y = gy(nu);
    minX = Math.min(minX, X);
    maxX = Math.max(maxX, X);
    minY = Math.min(minY, Y);
    maxY = Math.max(maxY, Y);
  }
  const mL = 26, // room for the "chief" label
    mR = 66, // room for the "perigee" label
    mT = 22,
    mB = g.perigee_window ? 34 : 24; // room for the caption
  const availW = W - mL - mR,
    availH = H - mT - mB;
  const k = Math.min(availW / (maxX - minX), availH / (maxY - minY));
  const figW = k * (maxX - minX),
    figH = k * (maxY - minY);
  const ox = mL + (availW - figW) / 2 - k * minX;
  const oy = mT + (availH - figH) / 2 + k * maxY;
  const px = (nu: number) => ox + k * gx(nu);
  const py = (nu: number) => oy - k * gy(nu);
  const figMidX = mL + availW / 2;

  // Orbit ellipse (sampled, focus on the chief).
  let d = "";
  for (let i = 0; i <= 180; i++) {
    const nu = (i / 180) * 2 * Math.PI - Math.PI;
    d += `${i === 0 ? "M" : "L"} ${px(nu).toFixed(2)} ${py(nu).toFixed(2)} `;
  }
  s.appendChild(el("path", { d: d + "Z", class: "orbit" }));
  // Perigee-window wedge (piecewise only).
  if (g.perigee_window) {
    const [lo, hi] = g.perigee_window;
    let wd = `M ${ox.toFixed(2)} ${oy.toFixed(2)} L ${px(lo).toFixed(2)} ${py(lo).toFixed(2)} `;
    const steps = 24;
    for (let i = 1; i <= steps; i++) {
      const nu = lo + ((hi - lo) * i) / steps;
      wd += `L ${px(nu).toFixed(2)} ${py(nu).toFixed(2)} `;
    }
    s.appendChild(el("path", { d: wd + "Z", class: "perigee-window" }));
  }
  // Chief (focus) + perigee marker, both labeled.
  s.appendChild(el("circle", { cx: ox, cy: oy, r: 4, class: "chief" }));
  s.appendChild(
    el("text", { x: ox - 8, y: oy + 4, class: "axis-label", "text-anchor": "end" }, "chief"),
  );
  s.appendChild(el("circle", { cx: px(0), cy: py(0), r: 3.5, class: "perigee" }));
  s.appendChild(
    el("text", { x: px(0) + 9, y: py(0) + 4, class: "axis-label" }, "perigee"),
  );
  // Maneuver markers — labels nudged radially outward so near-coincident burns
  // (and the markers themselves) don't collide.
  g.maneuver_nu.forEach((nu, j) => {
    const mx = px(nu),
      my = py(nu);
    s.appendChild(el("circle", { cx: mx, cy: my, r: 5, class: "mnvr-marker" }));
    s.appendChild(
      el(
        "text",
        {
          x: mx + Math.cos(nu) * 15,
          y: my - Math.sin(nu) * 15 + 4,
          class: "mnvr-label",
          "text-anchor": "middle",
        },
        String(j + 1),
      ),
    );
  });
  // Caption for the shaded wedge (piecewise perigee burn window).
  if (g.perigee_window) {
    s.appendChild(
      el(
        "text",
        { x: figMidX, y: H - 8, class: "axis-label", "text-anchor": "middle" },
        "shaded · perigee burn window",
      ),
    );
  }
  return s;
}

function panel(title: string, body: Node): HTMLElement {
  const p = document.createElement("section");
  p.className = "panel";
  const h = document.createElement("h2");
  h.textContent = title;
  p.append(h, body);
  return p;
}

export function render(outcome: SolveOutcome): void {
  const out = document.querySelector<HTMLElement>("#output")!;
  out.replaceChildren();
  // Reflect the solve outcome on the console status lamp (CSS-driven).
  document.body.classList.toggle("fault", outcome.status === "err");
  if (outcome.status === "err") {
    const err = document.createElement("div");
    err.className = `error ${outcome.error.kind}`;
    err.textContent = `${outcome.error.kind}: ${outcome.error.message}`;
    out.appendChild(err);
    return;
  }
  const r = outcome.value;
  const t_i = r.maneuvers.length ? Math.min(...r.maneuvers.map((m) => m.t)) : 0;
  const t_f = r.maneuvers.length ? Math.max(...r.maneuvers.map((m) => m.t)) : 1;
  out.appendChild(kpis(r));
  out.appendChild(panel("Δv timeline", timeline(r, t_i, t_f)));
  out.appendChild(panel("Δv components (R/T/N)", rtnBars(r)));
  out.appendChild(panel("Orbit geometry", orbit(r.geometry)));
}
