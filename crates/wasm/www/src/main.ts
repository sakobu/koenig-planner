import init, {
  solve,
  version,
  type SolveRequest,
  type CostSpec,
} from "koenig-damico-planner-wasm";
import { render } from "./render";
import { GOLDEN } from "./defaults";

function setNum(form: HTMLFormElement, name: string, value: number) {
  (form.elements.namedItem(name) as HTMLInputElement).value = String(value);
}

function getNum(form: HTMLFormElement, name: string): number {
  return Number((form.elements.namedItem(name) as HTMLInputElement).value);
}

function getOptNum(form: HTMLFormElement, name: string): number | undefined {
  const raw = (form.elements.namedItem(name) as HTMLInputElement).value.trim();
  return raw === "" ? undefined : Number(raw);
}

function writeForm(form: HTMLFormElement, r: SolveRequest) {
  setNum(form, "a", r.chief.a);
  setNum(form, "e", r.chief.e);
  setNum(form, "i", r.chief.i);
  setNum(form, "raan", r.chief.raan);
  setNum(form, "argp", r.chief.argp);
  setNum(form, "mean_anom", r.chief.mean_anom);
  setNum(form, "t_i", r.t_i);
  setNum(form, "t_f", r.t_f);
  setNum(form, "dt", r.dt);
  r.w_metres.forEach((w, k) => setNum(form, `w${k}`, w));
  (form.elements.namedItem("cost") as HTMLSelectElement).value = r.cost.type;
}

function readCost(form: HTMLFormElement): CostSpec {
  const type = (form.elements.namedItem("cost") as HTMLSelectElement).value;
  if (type === "norm2") return { type: "norm2" };
  if (type === "facemax") return { type: "facemax" };
  return {
    type: "piecewise",
    period: getOptNum(form, "period"),
    t_perigee0: getOptNum(form, "t_perigee0"),
  };
}

function readForm(form: HTMLFormElement): SolveRequest {
  // `params` (solver tuning) and `initial_times` (manual candidate-time
  // seeding) are deliberately omitted — the demo form exposes only the science
  // inputs and lets the core apply its defaults. Both are optional request
  // fields, so leaving them out is the documented default path.
  return {
    chief: {
      a: getNum(form, "a"),
      e: getNum(form, "e"),
      i: getNum(form, "i"),
      raan: getNum(form, "raan"),
      argp: getNum(form, "argp"),
      mean_anom: getNum(form, "mean_anom"),
    },
    t_i: getNum(form, "t_i"),
    t_f: getNum(form, "t_f"),
    dt: getNum(form, "dt"),
    w_metres: [
      getNum(form, "w0"),
      getNum(form, "w1"),
      getNum(form, "w2"),
      getNum(form, "w3"),
      getNum(form, "w4"),
      getNum(form, "w5"),
    ],
    cost: readCost(form),
  };
}

// Mirror the selected cost model onto the form's `data-cost` attribute so CSS
// can reveal the piecewise-only inputs (.pw) for piecewise alone.
function syncCost(form: HTMLFormElement) {
  form.dataset.cost = (
    form.elements.namedItem("cost") as HTMLSelectElement
  ).value;
}

function debounce<T extends (...a: never[]) => void>(fn: T, ms: number): T {
  let id: number | undefined;
  return ((...a: never[]) => {
    if (id !== undefined) clearTimeout(id);
    id = setTimeout(() => fn(...a), ms) as unknown as number;
  }) as T;
}

async function main() {
  await init();
  document.querySelector<HTMLElement>("#version")!.textContent =
    `core v${version()}`;
  const form = document.querySelector<HTMLFormElement>("#planner")!;
  writeForm(form, GOLDEN);
  const run = () => render(solve(readForm(form)));
  const debouncedRun = debounce(run, 150);
  form.addEventListener("input", () => {
    syncCost(form); // instant — toggle the piecewise-only inputs
    debouncedRun(); // debounced — re-solve
  });
  syncCost(form); // initial visibility (GOLDEN defaults to piecewise)
  run(); // initial render (instant)
}

main();
