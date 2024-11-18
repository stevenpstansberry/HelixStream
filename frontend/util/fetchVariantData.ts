/**
 * @file fetchVariantData.ts
 * @description This file contains the function to fetch variant data from an AWS S3 bucket.
 * @module fetchVariantData
 */

import AWS from 'aws-sdk';

AWS.config.update({
  region: process.env.NEXT_PUBLIC_AWS_REGION,
  accessKeyId: process.env.NEXT_PUBLIC_AWS_ACCESS_KEY_ID,
  secretAccessKey: process.env.NEXT_PUBLIC_AWS_SECRET_ACCESS_KEY,
});

const s3 = new AWS.S3();

/**
 * Fetches variant data from an AWS S3 bucket.
 * @param {string} variantName - The name of the variant to fetch data for.
 * @returns {Promise<any>} A promise that resolves to the variant data.
 */
export async function fetchVariantData(variantName: string): Promise<any> {
  const capitalizedVariantName = variantName.charAt(0).toUpperCase() + variantName.slice(1).toLowerCase();
  const params = {
    Bucket: 'bioinformantics-sequence-data',
    Key: `${capitalizedVariantName}/${variantName.toLowerCase()}.json`,
  };

  console.log(`Fetching from S3 bucket: ${params.Bucket}, file key: ${params.Key}`);

  return new Promise((resolve, reject) => {
    s3.getObject(params, (err, data) => {
      if (err) {
        reject(err);
      } else {
        if (data.Body) {
          resolve(JSON.parse(data.Body.toString('utf-8')));
        } else {
          reject(new Error('Data body is undefined'));
        }
      }
    });
  });
}