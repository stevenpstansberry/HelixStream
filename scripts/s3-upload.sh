#!/bin/bash

# Variables
BUCKET_NAME="your-s3-bucket-name"
LOCAL_DIR="../data/analyzed-sequences"
S3_DIR="s3://bioinformantics-sequence-data"

# Check if the local directory exists
if [ ! -d "$LOCAL_DIR" ]; then
  echo "Local directory $LOCAL_DIR does not exist."
  exit 1
fi

# Upload the directories to S3
echo "Uploading $LOCAL_DIR to $S3_DIR"
aws s3 cp "$LOCAL_DIR" "$S3_DIR" --recursive

# Check if the upload was successful
if [ $? -eq 0 ]; then
  echo "Upload successful."
else
  echo "Upload failed."
  exit 1
fi