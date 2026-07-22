import type { Metadata } from "next";
import { DesignSpaceApp } from "./design-space-app";

export const metadata: Metadata = {
  title: "SINEUP Delivery Atlas",
  description: "Disease-guided AAV and SINEUP delivery design space.",
};

export default function Home() {
  return <DesignSpaceApp />;
}
