"use client";

import { useEffect, useState } from "react";
import { Card } from "@/components/ui/card";
import { Info } from "lucide-react";
import AWS from 'aws-sdk';

export function InfoCard() {
  const [lastUpdated, setLastUpdated] = useState<string | null>(null);

  useEffect(() => {
    // Configure AWS SDK
    AWS.config.update({
      region: process.env.REACT_APP_AWS_REGION,
      accessKeyId: process.env.REACT_APP_AWS_ACCESS_KEY_ID,
      secretAccessKey: process.env.REACT_APP_AWS_SECRET_ACCESS_KEY,
    });

    const s3 = new AWS.S3();
    const params = {
      Bucket: 'bioinformantics-sequence-data',
      Key: 'lastupdated.json',
    };

    // Fetch the lastupdated.json file from S3
    s3.getObject(params, (err, data) => {
      if (err) {
        console.error('Error fetching lastupdated.json:', err);
        return;
      }

      if (data.Body) {
        const jsonData = JSON.parse(data.Body.toString('utf-8'));
        const lastUpdatedDate = new Date(jsonData.last_updated);

        // Format the date and time as a readable string
        const formattedDate = lastUpdatedDate.toLocaleString("en-US", {
          hour: "numeric",
          minute: "numeric",
          hour12: true,
          month: "long",
          day: "numeric",
          year: "numeric",
        });

        // Set the formatted date and time in the state
        setLastUpdated(formattedDate);
      } else {
        console.error('Error: data.Body is undefined');
      }
    });
  }, []);

  return (
    <Card className="p-6 mb-8">
      <div className="flex items-start gap-3">
        <Info className="h-5 w-5 text-blue-500 mt-1 flex-shrink-0" />
        <div>
          <p className="text-muted-foreground">
            This dashboard demonstrates a usecase for the HelixStream bioinformatics pipeline, designed to track the genetic evolution of SARS-CoV-2 variants over time. 
            Leveraging data from the NCBI database, the pipeline queues, ingests, and processes genome sequences across several prominent COVID-19 variants, 
            providing insights into mutation patterns and evolutionary trends. Key metrics highlight mutation rates and differences from reference genomes, 
            offering researchers valuable information on how these variants evolve and diverge.
            The pipleline will automatically update with new data as it becomes available and process it, ensuring that the most recent information is always accessible.
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