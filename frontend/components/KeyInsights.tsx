"use client";

import { Card } from "@/components/ui/card";
import { AlertCircle } from "lucide-react";

export function KeyInsights() {
    const data = [
        {
          Date: "2024-01-01",
          Percent_Identity: 98.03,
          Total_Mutations: 600,
          Mutation_Density_per_kb: 19.66,
          Jukes_Cantor_Distance: 0.0149,
          Substitutions: 440,
          Insertions: 148,
          Deletions: 12
        },
        {
          Date: "2024-01-15",
          Percent_Identity: 97.89,
          Total_Mutations: 620,
          Mutation_Density_per_kb: 20.12,
          Jukes_Cantor_Distance: 0.0156,
          Substitutions: 455,
          Insertions: 152,
          Deletions: 13
        },
        {
          Date: "2024-02-01",
          Percent_Identity: 97.76,
          Total_Mutations: 642,
          Mutation_Density_per_kb: 20.65,
          Jukes_Cantor_Distance: 0.0164,
          Substitutions: 470,
          Insertions: 157,
          Deletions: 15
        },
        {
          Date: "2024-02-15",
          Percent_Identity: 97.62,
          Total_Mutations: 665,
          Mutation_Density_per_kb: 21.23,
          Jukes_Cantor_Distance: 0.0173,
          Substitutions: 486,
          Insertions: 162,
          Deletions: 17
        },
        {
          Date: "2024-03-01",
          Percent_Identity: 97.45,
          Total_Mutations: 690,
          Mutation_Density_per_kb: 21.88,
          Jukes_Cantor_Distance: 0.0183,
          Substitutions: 504,
          Insertions: 168,
          Deletions: 18
        }
      ];
    
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
