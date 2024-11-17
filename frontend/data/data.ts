import { fetchVariantData } from '../util/fetchVariantData';

// Define the VariantData interface
interface VariantData {
  Date: string;
  Percent_Identity: number;
  Total_Mutations: number;
  Mutation_Density_per_kb: number;
  Jukes_Cantor_Distance: number;
  Substitutions: number;
  Insertions: number;
  Deletions: number;
}

// Omicron variant data
export let omicronData: VariantData[] = [];
export let omicronMutationTypeData: any[] = [];

fetchVariantData('omicron')
  .then((data: VariantData[]) => {
    omicronData = data;
    console.log('Omicron data fetched:', omicronData);

    // Set omicronMutationTypeData based on the latest record
    const latestOmicronRecord = omicronData[omicronData.length - 1];
    omicronMutationTypeData = [
      { name: 'Substitutions', value: latestOmicronRecord.Substitutions },
      { name: 'Insertions', value: latestOmicronRecord.Insertions },
      { name: 'Deletions', value: latestOmicronRecord.Deletions },
    ];
  })
  .catch((err: Error) => {
    console.error('Error fetching Omicron data:', err);
  });

// Alpha variant data
export let alphaData: VariantData[] = [];
export let alphaMutationTypeData: any[] = [];

fetchVariantData('alpha')
  .then((data: VariantData[]) => {
    alphaData = data;
    console.log('Alpha data fetched:', alphaData);

    // Set alphaMutationTypeData based on the latest record
    const latestAlphaRecord = alphaData[alphaData.length - 1];
    alphaMutationTypeData = [
      { name: 'Substitutions', value: latestAlphaRecord.Substitutions },
      { name: 'Insertions', value: latestAlphaRecord.Insertions },
      { name: 'Deletions', value: latestAlphaRecord.Deletions },
    ];
  })
  .catch((err: Error) => {
    console.error('Error fetching Alpha data:', err);
  });

// Beta variant data
export let betaData: VariantData[] = [];
export let betaMutationTypeData: any[] = [];

fetchVariantData('beta')
  .then((data: VariantData[]) => {
    betaData = data;
    console.log('Beta data fetched:', betaData);

    // Set betaMutationTypeData based on the latest record
    const latestBetaRecord = betaData[betaData.length - 1];
    betaMutationTypeData = [
      { name: 'Substitutions', value: latestBetaRecord.Substitutions },
      { name: 'Insertions', value: latestBetaRecord.Insertions },
      { name: 'Deletions', value: latestBetaRecord.Deletions },
    ];
  })
  .catch((err: Error) => {
    console.error('Error fetching Beta data:', err);
  });

// Gamma variant data
export let gammaData: VariantData[] = [];
export let gammaMutationTypeData: any[] = [];

fetchVariantData('gamma')
  .then((data: VariantData[]) => {
    gammaData = data;
    console.log('Gamma data fetched:', gammaData);

    // Set gammaMutationTypeData based on the latest record
    const latestGammaRecord = gammaData[gammaData.length - 1];
    gammaMutationTypeData = [
      { name: 'Substitutions', value: latestGammaRecord.Substitutions },
      { name: 'Insertions', value: latestGammaRecord.Insertions },
      { name: 'Deletions', value: latestGammaRecord.Deletions },
    ];
  })
  .catch((err: Error) => {
    console.error('Error fetching Gamma data:', err);
  });

// Delta variant data
export let deltaData: VariantData[] = [];
export let deltaMutationTypeData: any[] = [];

fetchVariantData('delta')
  .then((data: VariantData[]) => {
    deltaData = data;
    console.log('Delta data fetched:', deltaData);

    // Set deltaMutationTypeData based on the latest record
    const latestDeltaRecord = deltaData[deltaData.length - 1];
    deltaMutationTypeData = [
      { name: 'Substitutions', value: latestDeltaRecord.Substitutions },
      { name: 'Insertions', value: latestDeltaRecord.Insertions },
      { name: 'Deletions', value: latestDeltaRecord.Deletions },
    ];
  })
  .catch((err: Error) => {
    console.error('Error fetching Delta data:', err);
  });