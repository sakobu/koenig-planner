import type { ReactNode } from "react";

export function Panel({
  title,
  caption,
  children,
}: {
  title: string;
  caption?: ReactNode;
  children: ReactNode;
}) {
  return (
    <div className="panel">
      <h2>{title}</h2>
      {caption && <p className="panel-caption">{caption}</p>}
      {children}
    </div>
  );
}
