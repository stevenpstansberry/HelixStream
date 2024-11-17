"use client";

import { Card } from "@/components/ui/card";
import { AlertCircle } from "lucide-react";

type DataType = {
  Date: string;
  Percent_Identity: number;
  Total_Mutations: number;
  Mutation_Density_per_kb: number;
  Jukes_Cantor_Distance: number;
  Substitutions: number;
  Insertions: number;
  Deletions: number;
};

type KeyInsightsProps = {
  data: DataType[];
};

export function KeyInsights({ data }: KeyInsightsProps) {
  if (!data || data.length === 0) {
    return (
      <Card className="mt-8 p-6">
        <div className="flex items-start gap-3">
          <AlertCircle className="h-5 w-5 text-yellow-500 mt-1 flex-shrink-0" />
          <div>
            <h3 className="font-semibold mb-2">Key Insights</h3>
            <p className="text-muted-foreground">Loading data...</p>
          </div>
        </div>
      </Card>
    );
  }

  const latestData = data[data.length - 1];

  return (
    <Card className="mt-8 p-6">
      <div className="flex items-start gap-3">
        <AlertCircle className="h-5 w-5 text-yellow-500 mt-1 flex-shrink-0" />
        <div>
          <h3 className="font-semibold mb-2">Key Insights</h3>
          <p className="text-muted-foreground">
            Recent analysis shows a steady decline in sequence identity ({latestData.Percent_Identity.toFixed(2)}%), coupled with increased mutation density ({latestData.Mutation_Density_per_kb.toFixed(2)}/kb).
            Substitutions remain the dominant mutation type, accounting for {((latestData.Substitutions / latestData.Total_Mutations) * 100).toFixed(0)}% of all changes.
            The rising Jukes-Cantor distance ({latestData.Jukes_Cantor_Distance.toFixed(4)}) suggests continued evolutionary divergence, potentially affecting vaccine efficacy and viral characteristics.
          </p>
        </div>
      </div>
    </Card>
  );
}