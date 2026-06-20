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

function svg(w: number, h: number): SVGSVGElement {
  const s = el("svg", {
    viewBox: `0 0 ${w} ${h}`,
    width: "100%",
    preserveAspectRatio: "xMidYMid meet",
  }) as SVGSVGElement;
  return s;
}

function kpis(r: SolveResponse): HTMLElement {
  const box = document.createElement("div");
  box.className = "kpis";
  const cells: [string, string][] = [
    ["Σ Δv", `${r.total_dv.toFixed(4)} m/s`],
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
  const W = 640,
    H = 180,
    pad = 36;
  const s = svg(W, H);
  const mags = r.maneuvers.map((m) => Math.hypot(m.dv[0], m.dv[1], m.dv[2]));
  const maxMag = Math.max(1e-12, ...mags);
  const x = (t: number) =>
    pad + ((t - t_i) / Math.max(1e-9, t_f - t_i)) * (W - 2 * pad);
  const y = (mag: number) => H - pad - (mag / maxMag) * (H - 2 * pad);
  // Horizontal gridlines + magnitude ticks (0, ½·max, max).
  for (const frac of [0, 0.5, 1]) {
    const gy = y(frac * maxMag);
    s.appendChild(
      el("line", { x1: pad, y1: gy, x2: W - pad, y2: gy, class: "grid" }),
    );
    s.appendChild(
      el(
        "text",
        { x: pad - 6, y: gy + 3, class: "axis-label", "text-anchor": "end" },
        (frac * maxMag).toExponential(1),
      ),
    );
  }
  s.appendChild(
    el("line", {
      x1: pad,
      y1: H - pad,
      x2: W - pad,
      y2: H - pad,
      class: "axis",
    }),
  );
  s.appendChild(
    el(
      "text",
      {
        x: W - pad,
        y: H - pad + 16,
        class: "axis-label",
        "text-anchor": "end",
      },
      "t [s]",
    ),
  );
  s.appendChild(
    el(
      "text",
      { x: pad, y: H - pad + 16, class: "axis-label", "text-anchor": "start" },
      t_i.toFixed(0),
    ),
  );
  s.appendChild(
    el(
      "text",
      { x: pad - 6, y: 12, class: "axis-label", "text-anchor": "end" },
      "|Δv| [m/s]",
    ),
  );
  r.maneuvers.forEach((m, j) => {
    const mag = mags[j];
    s.appendChild(
      el("line", {
        x1: x(m.t),
        y1: H - pad,
        x2: x(m.t),
        y2: y(mag),
        class: "stem",
      }),
    );
    s.appendChild(
      el("circle", { cx: x(m.t), cy: y(mag), r: 4, class: "stem-dot" }),
    );
    s.appendChild(
      el(
        "text",
        {
          x: x(m.t),
          y: y(mag) - 8,
          class: "stem-label",
          "text-anchor": "middle",
        },
        mag.toExponential(1),
      ),
    );
  });
  return s;
}

function rtnBars(r: SolveResponse): SVGSVGElement {
  const W = 640,
    rowH = 28,
    pad = 90,
    H = r.maneuvers.length * rowH + 40;
  const s = svg(W, H);
  const maxComp = Math.max(
    1e-12,
    ...r.maneuvers.flatMap((m) => m.dv.map(Math.abs)),
  );
  const mid = (W + pad) / 2;
  const scale = (W - pad - 16 - (mid - pad)) / maxComp;
  // Zero axis spanning all rows.
  s.appendChild(
    el("line", {
      x1: mid,
      y1: 12,
      x2: mid,
      y2: 12 + r.maneuvers.length * rowH,
      class: "axis",
    }),
  );
  r.maneuvers.forEach((m, j) => {
    const yc = 12 + j * rowH + rowH / 2;
    s.appendChild(
      el("text", { x: 8, y: yc + 4, class: "row-label" }, `mnvr ${j + 1}`),
    );
    (["R", "T", "N"] as const).forEach((comp, k) => {
      const v = m.dv[k];
      const by = 12 + j * rowH + 3 + k * ((rowH - 6) / 3);
      const len = Math.abs(v) * scale;
      s.appendChild(
        el("rect", {
          x: v >= 0 ? mid : mid - len,
          y: by,
          width: len,
          height: (rowH - 6) / 3 - 2,
          fill: RTN_COLORS[comp],
        }),
      );
    });
  });
  // Component legend (R/T/N color key) along the bottom.
  (["R", "T", "N"] as const).forEach((comp, k) => {
    const lx = pad + k * 60;
    const ly = H - 14;
    s.appendChild(
      el("rect", {
        x: lx,
        y: ly - 8,
        width: 10,
        height: 10,
        fill: RTN_COLORS[comp],
      }),
    );
    s.appendChild(
      el("text", { x: lx + 15, y: ly, class: "legend-label" }, comp),
    );
  });
  return s;
}

function orbit(g: ChiefGeometry): SVGSVGElement {
  const W = 360,
    H = 360,
    cx = W / 2,
    cy = H / 2;
  const s = svg(W, H);
  const e = g.e,
    a = 1;
  const rOf = (nu: number) => (a * (1 - e * e)) / (1 + e * Math.cos(nu));
  // Fit: max radius is at apoapsis (ν = π).
  const rMax = rOf(Math.PI);
  const k = (Math.min(W, H) / 2 - 28) / rMax;
  const px = (nu: number) => cx + k * rOf(nu) * Math.cos(nu);
  const py = (nu: number) => cy - k * rOf(nu) * Math.sin(nu);
  // Orbit ellipse (sampled, focus-centered on the chief).
  let d = "";
  for (let i = 0; i <= 180; i++) {
    const nu = (i / 180) * 2 * Math.PI - Math.PI;
    d += `${i === 0 ? "M" : "L"} ${px(nu).toFixed(2)} ${py(nu).toFixed(2)} `;
  }
  s.appendChild(el("path", { d: d + "Z", class: "orbit" }));
  // Perigee-window wedge (piecewise only).
  if (g.perigee_window) {
    const [lo, hi] = g.perigee_window;
    let wd = `M ${cx} ${cy} L ${px(lo).toFixed(2)} ${py(lo).toFixed(2)} `;
    const steps = 24;
    for (let i = 1; i <= steps; i++) {
      const nu = lo + ((hi - lo) * i) / steps;
      wd += `L ${px(nu).toFixed(2)} ${py(nu).toFixed(2)} `;
    }
    s.appendChild(el("path", { d: wd + "Z", class: "perigee-window" }));
  }
  // Chief (focus) + perigee marker.
  s.appendChild(el("circle", { cx, cy, r: 4, class: "chief" }));
  s.appendChild(el("circle", { cx: px(0), cy: py(0), r: 3, class: "perigee" }));
  s.appendChild(
    el("text", { x: px(0) + 6, y: py(0) + 4, class: "axis-label" }, "perigee"),
  );
  // Maneuver markers.
  g.maneuver_nu.forEach((nu, j) => {
    s.appendChild(
      el("circle", { cx: px(nu), cy: py(nu), r: 5, class: "mnvr-marker" }),
    );
    s.appendChild(
      el(
        "text",
        { x: px(nu) + 7, y: py(nu) + 4, class: "mnvr-label" },
        String(j + 1),
      ),
    );
  });
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
