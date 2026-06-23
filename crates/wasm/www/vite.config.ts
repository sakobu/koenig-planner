import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import wasm from "vite-plugin-wasm";

export default defineConfig({
  base: "./",
  plugins: [react(), wasm()],
  server: {
    fs: {
      // The wasm package is a linked dependency (`file:../pkg`) whose real path
      // (crates/wasm/pkg) lives outside this demo root; allow the dev server to
      // serve it. (Production `vite build` bundles it, so this affects dev only.)
      allow: [".."],
    },
  },
});
