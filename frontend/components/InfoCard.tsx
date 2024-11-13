"use client";

import { useEffect, useState } from "react";
import { Card } from "@/components/ui/card";
import { Info } from "lucide-react";

export function InfoCard() {
  const [lastUpdated, setLastUpdated] = useState<string | null>(null);

  //TODO retrieve last-updated.json from bucket
  useEffect(() => {
    // Get the current date
    const currentDate = new Date();
    
    // Set the time to 3:00 PM
    currentDate.setHours(15, 0, 0, 0); // 15:00 in 24-hour format

    // Format the date and time as a readable string
    const formattedDate = currentDate.toLocaleString("en-US", {
      hour: "numeric",
      minute: "numeric",
      hour12: true,
      month: "long",
      day: "numeric",
      year: "numeric",
    });

    // Set the formatted date and time in the state
    setLastUpdated(formattedDate);
  }, []);

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
          {lastUpdated && (
            <p className="text-sm text-gray-500 mt-4">
              Last Updated: {lastUpdated}
            </p>
          )}
        </div>
      </div>
    </Card>
  );
}
