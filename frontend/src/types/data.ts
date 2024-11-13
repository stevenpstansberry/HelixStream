export interface SequenceData {
  Date: string;
  Sequence_ID: string;
  Percent_Identity: number;
  Substitutions: number;
  Insertions: number;
  Deletions: number;
  Total_Mutations: number;
  Mutation_Density_per_kb: number;
  Jukes_Cantor_Distance: number;
}