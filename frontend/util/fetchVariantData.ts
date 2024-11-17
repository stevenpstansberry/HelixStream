import AWS from 'aws-sdk';

AWS.config.update({
  region: process.env.NEXT_PUBLIC_AWS_REGION,
  accessKeyId: process.env.NEXT_PUBLIC_AWS_ACCESS_KEY_ID,
  secretAccessKey: process.env.NEXT_PUBLIC_AWS_SECRET_ACCESS_KEY,
});

const s3 = new AWS.S3();

export async function fetchVariantData(variantName: string): Promise<any> {
  const params = {
    Bucket: 'bioinformantics-sequence-data',
    Key: `${variantName.toLowerCase()}/${variantName.toLowerCase()}.json`,
  };

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