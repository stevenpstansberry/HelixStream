"use client";

import { useTheme } from "next-themes";
import { Dna, Activity, GitBranch, Microscope, Info, AlertCircle, Github, Linkedin } from "lucide-react"; // import icons
import { Header } from "@/components/Header";
import { InfoCard } from "@/components/InfoCard";
import { StatisticCard } from "@/components/StatisticCard";
import { ChartTabs } from "@/components/ChartTabs";
import { KeyInsights } from "@/components/KeyInsights";

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
          {/* Header */}
          <Header/>

          {/* Info Card */}
          <InfoCard/>
          
          {/* Statistic Cards */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 mb-8">
            <StatisticCard
              icon={<Activity className="h-5 w-5 text-blue-500" />}
              title="Sequence Identity"
              tooltip="Shows similarity to reference genome. Lower values indicate greater divergence from the original SARS-CoV-2 strain."
              value={`${latestData.Percent_Identity.toFixed(2)}%`}
              subtitle="vs. original SARS-CoV-2 Genome"
            />
            <StatisticCard
              icon={<GitBranch className="h-5 w-5 text-purple-500" />}
              title="Mutation Density"
              tooltip="Average mutations per kilobase, indicating the concentration of genetic changes in the viral genome."
              value={`${latestData.Mutation_Density_per_kb.toFixed(2)}/kb`}
              subtitle="Changes per Kilobase"
            />
            <StatisticCard
              icon={<Microscope className="h-5 w-5 text-pink-500" />}
              title="Total Mutations"
              tooltip="Aggregate count of all detected mutations, including substitutions, insertions, and deletions."
              value={`${latestData.Total_Mutations}`}
              subtitle="Combined Changes"
            />
            <StatisticCard
              icon={<Dna className="h-5 w-5 text-green-500" />}
              title="JC Distance"
              tooltip="Jukes-Cantor evolutionary distance, measuring genetic divergence while accounting for multiple mutations at the same site."
              value={`${latestData.Jukes_Cantor_Distance.toFixed(4)}`}
              subtitle="Evolutionary Distance"
            />
          </div>

          {/* Chart Tabs */}
          <ChartTabs/>

          {/* Key Insights */}
          <KeyInsights />
        </div>
        {/* Footer with social icons */}
        <footer className={`flex justify-center items-center py-4 border-t ${theme === 'dark' ? 'border-gray-600' : 'border-gray-300'}`}>
          <a href="https://github.com/stevenpstansberry" target="_blank" rel="noopener noreferrer" className="mx-4 group">
            <Github className="h-6 w-6 text-gray-600 group-hover:text-purple-500 group-hover:drop-shadow-[0_0_5px_rgba(128,0,128,0.75)] transition-all duration-300" />
          </a>
          <a href="https://linkedin.com/in/stevenpstansberry" target="_blank" rel="noopener noreferrer" className="mx-4 group">
            <Linkedin className="h-6 w-6 stroke-gray-600 fill-none group-hover:fill-purple-500 group-hover:drop-shadow-[0_0_5px_rgba(128,0,128,0.75)] transition-all duration-300" />
          </a>
        </footer>
      </div>
    </TooltipProvider>  
  );
}