"use client";

import { Card } from "@/components/ui/card";
import { Info } from "lucide-react";

export function InfoCard() {
  return (
    <Card className="p-6 mb-8">
      <div className="flex items-start gap-3">
        <Info className="h-5 w-5 text-blue-500 mt-1 flex-shrink-0" />
        <div>
          <p className="text-muted-foreground">
            This dashboard visualizes bioinformatics data processed through our custom pipeline designed to track the genetic evolution of SARS-CoV-2. 
            Focusing on mutations in the Omicron variant, we compare sequence changes over time to identify key differences and evolutionary trends. 
            The metrics displayed help researchers monitor how Omicron diverges from reference genomes, providing insights into mutation rates and potential functional impacts.
          </p>
        </div>
      </div>
    </Card>
  );
}
