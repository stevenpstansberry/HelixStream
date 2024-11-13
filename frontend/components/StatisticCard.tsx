"use client";

import { Card } from "@/components/ui/card";
import { Tooltip, TooltipContent, TooltipTrigger } from "@/components/ui/tooltip";
import { ReactNode } from "react";

interface StatisticCardProps {
  icon: ReactNode;
  title: string;
  tooltip: string;
  value: string;
  subtitle: string;
}

export function StatisticCard({ icon, title, tooltip, value, subtitle }: StatisticCardProps) {
  return (
    <Card className="p-6 hover:shadow-lg transition-shadow duration-300">
      <Tooltip>
        <TooltipTrigger asChild>
          <div className="flex items-center gap-3 cursor-help">
            {icon}
            <h3 className="text-sm font-medium text-muted-foreground">{title}</h3>
          </div>
        </TooltipTrigger>
        <TooltipContent side="top" sideOffset={4}>
          <p className="max-w-xs">{tooltip}</p>
        </TooltipContent>
      </Tooltip>
      <p className="text-2xl font-bold mt-2">{value}</p>
      <p className="text-sm text-muted-foreground mt-1">{subtitle}</p>
    </Card>
  );
}
