#!/bin/bash

# Variables
BUCKET_NAME="your-s3-bucket-name"
LOCAL_DIR="../data/analyzed-sequences"
S3_DIR="s3://bioinformantics-sequence-data"
LAST_UPDATED_FILE="../data/lastupdated.json"

# Check if the local directory exists
if [ ! -d "$LOCAL_DIR" ]; then
  echo "Local directory $LOCAL_DIR does not exist."
  exit 1
fi

# Create or overwrite the lastupdated.json file with the current timestamp
echo "Creating or overwriting lastupdated.json with the current timestamp..."
current_time=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
echo "{\"last_updated\": \"$current_time\"}" > "$LAST_UPDATED_FILE"

# Upload the directories to S3
echo "Uploading $LOCAL_DIR to $S3_DIR"
aws s3 cp "$LOCAL_DIR" "$S3_DIR" --recursive

# Check if the upload was successful
if [ $? -eq 0 ]; then
  echo "Upload successful."

  # Upload the lastupdated.json file to S3
  echo "Uploading lastupdated.json to $S3_DIR"
  aws s3 cp "$LAST_UPDATED_FILE" "$S3_DIR/lastupdated.json"

  # Check if the lastupdated.json upload was successful
  if [ $? -eq 0 ]; then
    echo "lastupdated.json upload successful."
  else
    echo "Failed to upload lastupdated.json."
    exit 1
  fi

  # Clean up files in each directory
  find "$LOCAL_DIR" -type f -exec rm -f {} \;
  echo "Files in $LOCAL_DIR have been removed."
else
  echo "Upload failed."
  exit 1
fi