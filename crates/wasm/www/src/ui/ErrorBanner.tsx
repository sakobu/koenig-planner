/** The single typed-error surface (wasm-init failure, and the solver /
 *  bad_request / internal API errors). `kind` drives BOTH the CSS variant
 *  (amber `bad_request` warning vs red fault) and the `kind: message` prefix;
 *  omit it for a bare message (init failure). `variant` adds the `internal`
 *  label or the sticky `overlay` treatment used above a still-mounted readout.
 *  `role="alert"` announces the fault to assistive tech — the one place that
 *  concern now lives. */
export function ErrorBanner({
  kind,
  message,
  variant,
}: {
  kind?: string;
  message: string;
  variant?: "overlay" | "internal";
}) {
  const className = ["error", kind, variant].filter(Boolean).join(" ");
  return (
    <div className={className} role="alert">
      {kind ? `${kind}: ${message}` : message}
    </div>
  );
}
