import React from 'react';
import { SequenceData } from '../types/data';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer
} from 'recharts';
import { format } from 'date-fns';

interface MutationBarChartProps {
  data: SequenceData[];
}

export const MutationBarChart: React.FC<MutationBarChartProps> = ({ data }) => {
  const formattedData = data.map(item => ({
    date: format(new Date(item.Date), 'MMM d, yyyy'),
    Substitutions: item.Substitutions,
    Insertions: item.Insertions,
    Deletions: item.Deletions
  }));

  return (
    <div className="bg-white rounded-lg shadow-md p-4 h-[400px]">
      <h2 className="text-lg font-semibold mb-4">Mutation Composition</h2>
      <ResponsiveContainer width="100%" height="100%">
        <BarChart data={formattedData} margin={{ top: 20, right: 30, left: 20, bottom: 5 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="date" />
          <YAxis />
          <Tooltip />
          <Legend />
          <Bar dataKey="Substitutions" stackId="a" fill="#3B82F6" />
          <Bar dataKey="Insertions" stackId="a" fill="#10B981" />
          <Bar dataKey="Deletions" stackId="a" fill="#EF4444" />
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
};