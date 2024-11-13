"use client";

import { Card } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, BarChart, Bar, PieChart, Pie, Cell } from "recharts";
import { ThemeToggle } from "@/components/theme-toggle";
import { useTheme } from "next-themes";
import { Dna, Activity, GitBranch, Microscope, Info, AlertCircle } from "lucide-react";
import { Button } from "@/components/ui/button";
import {
  Tooltip as UITooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";

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

const mutationTypeData = [
  { name: 'Substitutions', value: 504 },
  { name: 'Insertions', value: 168 },
  { name: 'Deletions', value: 18 }
];

const COLORS = ['hsl(var(--chart-2))', 'hsl(var(--chart-3))', 'hsl(var(--chart-4))'];

export default function Dashboard() {
  const { theme } = useTheme();
  const chartTheme = theme === 'dark' ? {
    backgroundColor: 'transparent',
    textColor: '#ffffff',
    gridColor: '#333333'
  } : {
    backgroundColor: 'transparent',
    textColor: '#000000',
    gridColor: '#e5e5e5'
  };
  
  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleDateString();
  };

  const latestData = data[data.length - 1];

  return (
    <TooltipProvider>
      <div className="min-h-screen bg-background transition-colors duration-300">
        <div className="mx-auto max-w-7xl p-6">
          <div className="flex justify-between items-center mb-6">
            <div className="flex items-center gap-3">
              <Dna className="h-8 w-8 text-primary" />
              <div>
                <h1 className="text-4xl font-bold">Bioinformatics Dashboard</h1>
                <h2 className="text-xl text-muted-foreground mt-1">COVID-19 Variant Mutation Analysis</h2>
              </div>
            </div>
            <ThemeToggle />
          </div>

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
          
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 mb-8">
            <TooltipProvider>
              <Card className="p-6 hover:shadow-lg transition-shadow duration-300">
                <UITooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center gap-3 cursor-help">
                      <Activity className="h-5 w-5 text-blue-500" />
                      <h3 className="text-sm font-medium text-muted-foreground">Sequence Identity</h3>
                    </div>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p className="max-w-xs">Shows similarity to reference genome. Lower values indicate greater divergence from the original SARS-CoV-2 strain.</p>
                  </TooltipContent>
                </UITooltip>
                <p className="text-2xl font-bold mt-2">{latestData.Percent_Identity.toFixed(2)}%</p>
                <p className="text-sm text-muted-foreground mt-1">vs. Reference Genome</p>
              </Card>

              <Card className="p-6 hover:shadow-lg transition-shadow duration-300">
                <UITooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center gap-3 cursor-help">
                      <GitBranch className="h-5 w-5 text-purple-500" />
                      <h3 className="text-sm font-medium text-muted-foreground">Mutation Density</h3>
                    </div>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p className="max-w-xs">Average mutations per kilobase, indicating the concentration of genetic changes in the viral genome.</p>
                  </TooltipContent>
                </UITooltip>
                <p className="text-2xl font-bold mt-2">{latestData.Mutation_Density_per_kb.toFixed(2)}/kb</p>
                <p className="text-sm text-muted-foreground mt-1">Changes per Kilobase</p>
              </Card>

              <Card className="p-6 hover:shadow-lg transition-shadow duration-300">
                <UITooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center gap-3 cursor-help">
                      <Microscope className="h-5 w-5 text-pink-500" />
                      <h3 className="text-sm font-medium text-muted-foreground">Total Mutations</h3>
                    </div>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p className="max-w-xs">Aggregate count of all detected mutations, including substitutions, insertions, and deletions.</p>
                  </TooltipContent>
                </UITooltip>
                <p className="text-2xl font-bold mt-2">{latestData.Total_Mutations}</p>
                <p className="text-sm text-muted-foreground mt-1">Combined Changes</p>
              </Card>

              <Card className="p-6 hover:shadow-lg transition-shadow duration-300">
                <UITooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center gap-3 cursor-help">
                      <Dna className="h-5 w-5 text-green-500" />
                      <h3 className="text-sm font-medium text-muted-foreground">JC Distance</h3>
                    </div>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p className="max-w-xs">Jukes-Cantor evolutionary distance, measuring genetic divergence while accounting for multiple mutations at the same site.</p>
                  </TooltipContent>
                </UITooltip>
                <p className="text-2xl font-bold mt-2">{latestData.Jukes_Cantor_Distance.toFixed(4)}</p>
                <p className="text-sm text-muted-foreground mt-1">Evolutionary Distance</p>
              </Card>
            </TooltipProvider>
          </div>

          <Tabs defaultValue="identity" className="space-y-4">
            <TabsList className="w-full justify-start">
              <TabsTrigger value="identity">Sequence Identity</TabsTrigger>
              <TabsTrigger value="mutations">Mutation Types</TabsTrigger>
              <TabsTrigger value="density">Mutation Density</TabsTrigger>
              <TabsTrigger value="distance">Evolutionary Distance</TabsTrigger>
            </TabsList>

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
                        tickFormatter={formatDate}
                        stroke={chartTheme.textColor}
                      />
                      <YAxis 
                        domain={[95, 100]}
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

            <TabsContent value="mutations">
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
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
                          tickFormatter={formatDate}
                          stroke={chartTheme.textColor}
                        />
                        <YAxis stroke={chartTheme.textColor} />
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
                        tickFormatter={formatDate}
                        stroke={chartTheme.textColor}
                      />
                      <YAxis stroke={chartTheme.textColor} />
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
                        tickFormatter={formatDate}
                        stroke={chartTheme.textColor}
                      />
                      <YAxis stroke={chartTheme.textColor} />
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
        </div>
      </div>
    </TooltipProvider>  
  );
}