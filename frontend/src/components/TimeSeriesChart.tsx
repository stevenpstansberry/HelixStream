import React from 'react';
import { SequenceData } from '../types/data';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer
} from 'recharts';
import { format } from 'date-fns';

interface TimeSeriesChartProps {
  data: SequenceData[];
  metric: keyof SequenceData;
  title: string;
  color: string;
}

export const TimeSeriesChart: React.FC<TimeSeriesChartProps> = ({
  data,
  metric,
  title,
  color
}) => {
  const formattedData = data.map(item => ({
    ...item,
    Date: format(new Date(item.Date), 'MMM d, yyyy')
  }));

  return (
    <div className="bg-white rounded-lg shadow-md p-4 h-[400px]">
      <h2 className="text-lg font-semibold mb-4">{title}</h2>
      <ResponsiveContainer width="100%" height="100%">
        <LineChart data={formattedData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="Date" />
          <YAxis />
          <Tooltip />
          <Legend />
          <Line
            type="monotone"
            dataKey={metric}
            stroke={color}
            strokeWidth={2}
            dot={{ r: 4 }}
            activeDot={{ r: 8 }}
          />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
};