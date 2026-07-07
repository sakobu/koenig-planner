// A full-width divider that titles a group of panels in the readout deck, so the
// family boundaries (geometry / plan / certificate / trade) read as sections
// rather than a flat wall of identical instrument cards.
export function SectionBand({ label, hint }: { label: string; hint?: string }) {
  return (
    <div className="section-band">
      <span className="band-label">{label}</span>
      {hint && <span className="band-hint">{hint}</span>}
    </div>
  );
}
