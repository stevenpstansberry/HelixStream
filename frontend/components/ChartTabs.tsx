"use client";

import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Card } from "@/components/ui/card";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer,
  BarChart, Bar, PieChart, Pie, Cell,
} from "recharts";
import { Button } from "@/components/ui/button";
import { Info } from "lucide-react";
import {
  Tooltip as UITooltip,
  TooltipContent,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import { useTheme } from "next-themes";

const mutationTypeData = [
  { name: 'Substitutions', value: 504 },
  { name: 'Insertions', value: 168 },
  { name: 'Deletions', value: 18 }
];

const COLORS = ['hsl(var(--chart-2))', 'hsl(var(--chart-3))', 'hsl(var(--chart-4))'];

const formatDate = (dateString: string) => {
  return new Date(dateString).toLocaleDateString();
};

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

const data: DataType[] = [
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

export function ChartTabs() {
  const { theme } = useTheme();
  const chartTheme = theme === 'dark' ? {
    textColor: '#ffffff',
    gridColor: '#333333'
  } : {
    backgroundColor: 'transparent',
    textColor: '#000000',
    gridColor: '#e5e5e5'
  };

  // Dynamically calculate Y-axis domains based on data
  const getYAxisDomain = (dataKey: string, buffer: number = 1) => {
    const values = data.map((d) => d[dataKey as keyof DataType]);
    const min = Math.min(...values.map(Number));
    const max = Math.max(...values.map(Number));
    return [Math.floor(min) - buffer, Math.ceil(max) + buffer];
  };

  return (
    <Tabs defaultValue="identity" className="space-y-4">
      <TabsList className="w-full justify-start">
        <TabsTrigger value="identity">Sequence Identity</TabsTrigger>
        <TabsTrigger value="mutations">Mutation Types</TabsTrigger>
        <TabsTrigger value="density">Mutation Density</TabsTrigger>
        <TabsTrigger value="distance">Evolutionary Distance</TabsTrigger>
      </TabsList>

      {/* Sequence Identity Tab */}
      <TabsContent value="identity">
        <Card className="p-6">
          <div className="flex items-center justify-between mb-4">
            <h2 className="text-xl font-semibold">Sequence Identity Analysis</h2>
            <UITooltip>
              <TooltipTrigger asChild>
                <Button variant="ghost" size="icon">
                  <Info className="h-4 w-4" />
                </Button>
              </TooltipTrigger>
              <TooltipContent>
                <p className="max-w-xs">Tracks similarity between Omicron sequences and the reference genome over time. Declining values suggest viral evolution.</p>
              </TooltipContent>
            </UITooltip>
          </div>
          <div className="h-[400px]">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data}>
                <CartesianGrid strokeDasharray="3 3" stroke={chartTheme.gridColor} />
                <XAxis
                    dataKey="Date"
                    label={{ value: 'Date', position: 'insideBottomRight', offset: 10, dy: 10, dx: -525 }}
                    tickFormatter={formatDate}
                    stroke={chartTheme.textColor}
                    />
                    <YAxis
                    label={{
                        value: 'Sequence Identity (%)',
                        angle: -90,
                        position: 'insideBottomLeft',
                        dx: 20,
                        dy: -75,
                    }}
                    domain={getYAxisDomain('Percent_Identity')}
                    stroke={chartTheme.textColor}
                />
                <Tooltip
                  labelFormatter={formatDate}
                  contentStyle={{
                    backgroundColor: theme === 'dark' ? '#1f2937' : '#ffffff',
                    borderColor: theme === 'dark' ? '#374151' : '#e5e5e5',
                    color: chartTheme.textColor
                  }}
                />
                <Legend />
                <Line
                  type="monotone"
                  dataKey="Percent_Identity"
                  stroke="hsl(var(--chart-1))"
                  name="Sequence Identity (%)"
                  strokeWidth={2}
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </Card>
      </TabsContent>

      {/* Mutation Types Tab */}
      <TabsContent value="mutations">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Mutation Distribution Bar Chart */}
          <Card className="p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-xl font-semibold">Mutation Distribution</h2>
              <UITooltip>
                <TooltipTrigger asChild>
                  <Button variant="ghost" size="icon">
                    <Info className="h-4 w-4" />
                  </Button>
                </TooltipTrigger>
                <TooltipContent>
                  <p className="max-w-xs">Breakdown of different mutation types over time, showing the composition of genetic changes.</p>
                </TooltipContent>
              </UITooltip>
            </div>
            <div className="h-[400px]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={data}>
                  <CartesianGrid strokeDasharray="3 3" stroke={chartTheme.gridColor} />
                  <XAxis
                    dataKey="Date"
                    label={{ value: 'Date', position: 'insideBottom', offset: -5, dx: -25}}
                    tickFormatter={formatDate}
                    stroke={chartTheme.textColor}
                  />
                  <YAxis
                    label={{ value: 'Number of Mutations', angle: -90, position: 'insideLeft', dy: 75 }}
                    domain={[0, 'auto']}
                    stroke={chartTheme.textColor}
                  />
                  <Tooltip
                    labelFormatter={formatDate}
                    contentStyle={{
                      backgroundColor: theme === 'dark' ? '#1f2937' : '#ffffff',
                      borderColor: theme === 'dark' ? '#374151' : '#e5e5e5',
                      color: chartTheme.textColor
                    }}
                  />
                  <Legend />
                  <Bar dataKey="Substitutions" stackId="a" fill="hsl(var(--chart-2))" name="Substitutions" />
                  <Bar dataKey="Insertions" stackId="a" fill="hsl(var(--chart-3))" name="Insertions" />
                  <Bar dataKey="Deletions" stackId="a" fill="hsl(var(--chart-4))" name="Deletions" />
                </BarChart>
              </ResponsiveContainer>
            </div>
          </Card>

          {/* Current Mutation Composition Pie Chart */}
          <Card className="p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-xl font-semibold">Current Mutation Composition</h2>
              <UITooltip>
                <TooltipTrigger asChild>
                  <Button variant="ghost" size="icon">
                    <Info className="h-4 w-4" />
                  </Button>
                </TooltipTrigger>
                <TooltipContent>
                  <p className="max-w-xs">Current distribution of mutation types in the latest sequence analysis.</p>
                </TooltipContent>
              </UITooltip>
            </div>
            <div className="h-[400px]">
              <ResponsiveContainer width="100%" height="100%">
                <PieChart>
                  <Pie
                    data={mutationTypeData}
                    cx="50%"
                    cy="50%"
                    labelLine={false}
                    label={({ name, percent }) => `${name} (${(percent * 100).toFixed(0)}%)`}
                    outerRadius={150}
                    fill="#8884d8"
                    dataKey="value"
                  >
                    {mutationTypeData.map((entry, index) => (
                      <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                    ))}
                  </Pie>
                  <Tooltip />
                </PieChart>
              </ResponsiveContainer>
            </div>
          </Card>
        </div>
      </TabsContent>

      {/* Mutation Density Tab */}
      <TabsContent value="density">
        <Card className="p-6">
          <div className="flex items-center justify-between mb-4">
            <h2 className="text-xl font-semibold">Mutation Density Analysis</h2>
            <UITooltip>
              <TooltipTrigger asChild>
                <Button variant="ghost" size="icon">
                  <Info className="h-4 w-4" />
                </Button>
              </TooltipTrigger>
              <TooltipContent>
                <p className="max-w-xs">Tracks mutation density over time, indicating the rate of genetic changes per kilobase.</p>
              </TooltipContent>
            </UITooltip>
          </div>
          <div className="h-[400px]">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data}>
                <CartesianGrid strokeDasharray="3 3" stroke={chartTheme.gridColor} />
                <XAxis
                  dataKey="Date"
                  label={{ value: 'Date', position: 'insideBottom', offset: 0, dx: -10 }}
                  tickFormatter={formatDate}
                  stroke={chartTheme.textColor}
                />
                <YAxis
                  label={{ value: 'Mutation Density (per kb)', angle: -90, position: 'insideLeft', dy: 75 }}
                  domain={getYAxisDomain('Mutation_Density_per_kb')}
                  stroke={chartTheme.textColor}
                />
                <Tooltip
                  labelFormatter={formatDate}
                  contentStyle={{
                    backgroundColor: theme === 'dark' ? '#1f2937' : '#ffffff',
                    borderColor: theme === 'dark' ? '#374151' : '#e5e5e5',
                    color: chartTheme.textColor
                  }}
                />
                <Legend />
                <Line
                  type="monotone"
                  dataKey="Mutation_Density_per_kb"
                  stroke="hsl(var(--chart-3))"
                  name="Mutations per Kilobase"
                  strokeWidth={2}
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </Card>
      </TabsContent>

      {/* Evolutionary Distance Tab */}
      <TabsContent value="distance">
        <Card className="p-6">
          <div className="flex items-center justify-between mb-4">
            <h2 className="text-xl font-semibold">Evolutionary Distance Trends</h2>
            <UITooltip>
              <TooltipTrigger asChild>
                <Button variant="ghost" size="icon">
                  <Info className="h-4 w-4" />
                </Button>
              </TooltipTrigger>
              <TooltipContent>
                <p className="max-w-xs">Jukes-Cantor distance showing evolutionary divergence from the reference genome.</p>
              </TooltipContent>
            </UITooltip>
          </div>
          <div className="h-[400px]">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data}>
                <CartesianGrid strokeDasharray="3 3" stroke={chartTheme.gridColor} />
                <XAxis
                  dataKey="Date"
                  label={{ value: 'Date', position: 'insideBottom', offset: -5 }}
                  tickFormatter={formatDate}
                  stroke={chartTheme.textColor}
                />
                <YAxis
                  label={{ value: 'Jukes-Cantor Distance', angle: -90, position: 'insideLeft', dy: 75 }}
                  domain={[0, 'auto']}
                  stroke={chartTheme.textColor}
                />
                <Tooltip
                  labelFormatter={formatDate}
                  contentStyle={{
                    backgroundColor: theme === 'dark' ? '#1f2937' : '#ffffff',
                    borderColor: theme === 'dark' ? '#374151' : '#e5e5e5',
                    color: chartTheme.textColor
                  }}
                />
                <Legend />
                <Line
                  type="monotone"
                  dataKey="Jukes_Cantor_Distance"
                  stroke="hsl(var(--chart-4))"
                  name="Evolutionary Distance"
                  strokeWidth={2}
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </Card>
      </TabsContent>
    </Tabs>
  );
}
