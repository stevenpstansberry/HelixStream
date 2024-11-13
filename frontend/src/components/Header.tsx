import React from 'react';
import { SequenceData } from '../types/data';
import { DNA, Activity, GitBranch } from 'lucide-react';

interface HeaderProps {
  data: SequenceData[];
}

export const Header: React.FC<HeaderProps> = ({ data }) => {
  const averageIdentity = data.reduce((acc, curr) => acc + curr.Percent_Identity, 0) / data.length;
  const averageDensity = data.reduce((acc, curr) => acc + curr.Mutation_Density_per_kb, 0) / data.length;
  const latestData = data[data.length - 1];

  return (
    <header className="bg-gradient-to-r from-blue-600 to-indigo-700 text-white p-6 shadow-lg">
      <div className="max-w-7xl mx-auto">
        <div className="flex items-center justify-between mb-6">
          <div className="flex items-center space-x-3">
            <DNA className="h-8 w-8" />
            <h1 className="text-2xl font-bold">Bioinformatics Dashboard</h1>
          </div>
          <div className="text-sm">
            Last updated: {new Date(latestData.Date).toLocaleString()}
          </div>
        </div>
        
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <div className="bg-white/10 rounded-lg p-4 backdrop-blur-sm">
            <div className="flex items-center space-x-2">
              <Activity className="h-5 w-5" />
              <h3 className="text-sm font-medium">Avg. Percent Identity</h3>
            </div>
            <p className="text-2xl font-bold mt-2">{averageIdentity.toFixed(2)}%</p>
          </div>
          
          <div className="bg-white/10 rounded-lg p-4 backdrop-blur-sm">
            <div className="flex items-center space-x-2">
              <GitBranch className="h-5 w-5" />
              <h3 className="text-sm font-medium">Avg. Mutation Density</h3>
            </div>
            <p className="text-2xl font-bold mt-2">{averageDensity.toFixed(2)} per kb</p>
          </div>
          
          <div className="bg-white/10 rounded-lg p-4 backdrop-blur-sm">
            <div className="flex items-center space-x-2">
              <DNA className="h-5 w-5" />
              <h3 className="text-sm font-medium">Total Sequences</h3>
            </div>
            <p className="text-2xl font-bold mt-2">{data.length}</p>
          </div>
        </div>
      </div>
    </header>
  );
};