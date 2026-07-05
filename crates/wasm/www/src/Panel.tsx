import type { ReactNode } from "react";

export function Panel({ title, children }: { title: string; children: ReactNode }) {
  return (
    <div className="panel">
      <h2>{title}</h2>
      {children}
    </div>
  );
}
