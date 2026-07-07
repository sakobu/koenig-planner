import { Component, type ReactNode } from "react";

/** Last-resort catch for render faults anywhere below <App/>: shows the message
 *  in the native error style instead of letting React unmount to a blank page. */
export class ErrorBoundary extends Component<
  { children: ReactNode },
  { error: Error | null }
> {
  state = { error: null as Error | null };

  static getDerivedStateFromError(error: Error) {
    return { error };
  }

  render() {
    if (this.state.error) {
      return (
        <section id="output">
          <div className="error internal">
            {`render fault: ${this.state.error.message}`}
          </div>
        </section>
      );
    }
    return this.props.children;
  }
}
