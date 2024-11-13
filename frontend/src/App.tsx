import React from 'react';
import { Header } from './components/Header';
import { TimeSeriesChart } from './components/TimeSeriesChart';
import { MutationBarChart } from './components/MutationBarChart';
import { SequenceData } from './types/data';

const sampleData: SequenceData[] = [
  {
    "Date": "2022-09-29T00:00:00.000",
    "Sequence_ID": "human/FRA/HEGP-81-10-2332993394i_p2/2022###2022-09-29",
    "Percent_Identity": 99.4350660204,
    "Substitutions": 94,
    "Insertions": 87,
    "Deletions": 12,
    "Total_Mutations": 181,
    "Mutation_Density_per_kb": 6.0504763497,
    "Jukes_Cantor_Distance": 0.0031501035
  },
  {
    "Date": "2023-08-21T00:00:00.000",
    "Sequence_ID": "human/FRA/IPP17625i_p1/2023###2023-08-21",
    "Percent_Identity": 98.2884840381,
    "Substitutions": 124,
    "Insertions": 388,
    "Deletions": 12,
    "Total_Mutations": 512,
    "Mutation_Density_per_kb": 17.1151596189,
    "Jukes_Cantor_Distance": 0.0041582472
  },
  {
    "Date": "2023-12-12T00:00:00.000",
    "Sequence_ID": "human/FRA/045-0010-C-D_2401V120212/2023###2023-12-12",
    "Percent_Identity": 98.6528497409,
    "Substitutions": 131,
    "Insertions": 272,
    "Deletions": 12,
    "Total_Mutations": 403,
    "Mutation_Density_per_kb": 13.4715025907,
    "Jukes_Cantor_Distance": 0.0043936759
  },
  {
    "Date": "2023-12-29T00:00:00.000",
    "Sequence_ID": "human/FRA/035-0012-P-G_2401V120198/2023###2023-12-29",
    "Percent_Identity": 98.9837873976,
    "Substitutions": 141,
    "Insertions": 163,
    "Deletions": 12,
    "Total_Mutations": 304,
    "Mutation_Density_per_kb": 10.1621260237,
    "Jukes_Cantor_Distance": 0.0047301307
  },
  {
    "Date": "2024-01-05T00:00:00.000",
    "Sequence_ID": "human/FRA/045-0011-D-P_2401V120213/2024###2024-01-05",
    "Percent_Identity": 98.0344308875,
    "Substitutions": 440,
    "Insertions": 148,
    "Deletions": 12,
    "Total_Mutations": 588,
    "Mutation_Density_per_kb": 19.6556911249,
    "Jukes_Cantor_Distance": 0.0148604981
  }
];

function App() {
  return (
    <div className="min-h-screen bg-gray-50">
      <Header data={sampleData} />
      
      <main className="max-w-7xl mx-auto px-4 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <TimeSeriesChart
            data={sampleData}
            metric="Percent_Identity"
            title="Percent Identity Over Time"
            color="#3B82F6"
          />
          
          <TimeSeriesChart
            data={sampleData}
            metric="Jukes_Cantor_Distance"
            title="Jukes-Cantor Distance Over Time"
            color="#8B5CF6"
          />
          
          <TimeSeriesChart
            data={sampleData}
            metric="Mutation_Density_per_kb"
            title="Mutation Density per kb Over Time"
            color="#10B981"
          />
          
          <MutationBarChart data={sampleData} />
        </div>
      </main>
    </div>
  );
}

export default App;