"use client";

import { useState } from "react";
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

import { alphaData, betaData, gammaData, deltaData, omicronData, alphaMutationTypeData, betaMutationTypeData, gammaMutationTypeData, deltaMutationTypeData, omicronMutationTypeData } from "@/data/data";


const variantOptions = {
  Alpha: alphaData,
  Beta: betaData,
  Gamma: gammaData,
  Delta: deltaData,
  Omicron: omicronData
};

const mutationTypeVariantOptions = {
  Alpha: alphaMutationTypeData,
  Beta: betaMutationTypeData,
  Gamma: gammaMutationTypeData,
  Delta: deltaMutationTypeData,
  Omicron: omicronMutationTypeData
};

export default function Dashboard() {
  const { theme } = useTheme();
  const [selectedVariant, setSelectedVariant] = useState("Omicron");
  const [variantData, setVariantData] = useState(variantOptions["Omicron"]);
  const [mutationTypeData, setMutationTypeData] = useState(omicronMutationTypeData);
  const chartTheme = theme === 'dark' ? {
    backgroundColor: 'transparent',
    textColor: '#ffffff',
    gridColor: '#333333'
  } : {
    backgroundColor: 'transparent',
    textColor: '#000000',
    gridColor: '#e5e5e5'
  };

  const handleVariantChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    const selected = event.target.value;
    setSelectedVariant(selected);
    setVariantData(variantOptions[selected as keyof typeof variantOptions]);
    setMutationTypeData(mutationTypeVariantOptions[selected as keyof typeof mutationTypeVariantOptions]);
  };

  
  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleDateString();
  };

  const latestData = variantData[variantData.length - 1];

  return (
    <TooltipProvider>
      <div className="min-h-screen bg-background transition-colors duration-300">
        <div className="mx-auto max-w-7xl p-6">
          {/* Header */}
          <Header/>


          {/* Variant Selection Dropdown */}
          <div className="mb-6">
            <label htmlFor="variantSelect" className="block text-sm font-medium text-gray-700">
              Select Variant:
            </label>
            <select
              id="variantSelect"
              value={selectedVariant}
              onChange={handleVariantChange}
              className="mt-1 block w-full p-2 border border-gray-300 rounded-md shadow-sm focus:border-indigo-500 focus:ring-indigo-500 sm:text-sm"
            >
              {Object.keys(variantOptions).map(variant => (
                <option key={variant} value={variant}>
                  {variant}
                </option>
              ))}
            </select>
          </div>

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
          <ChartTabs data={variantData} mutationTypeData={mutationTypeData}/>

          {/* Key Insights */}
          <KeyInsights data={variantData}/>
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