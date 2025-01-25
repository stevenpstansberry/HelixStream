#!/bin/bash

# Colors for messages
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Variables
LOCAL_DIR="../data/analyzed-sequences"
S3_DIR="s3://bioinformantics-sequence-data"
LAST_UPDATED_FILE="../data/lastupdated.json"
LOG_FILE="../logs/upload_$(date +%Y%m%d%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Trap interrupt signals for cleanup
trap 'echo -e "${RED}Process interrupted. Exiting...${NC}"; exit 1' INT

# Check if the local directory exists
if [ ! -d "$LOCAL_DIR" ]; then
  echo -e "${RED}Local directory $LOCAL_DIR does not exist.${NC}"
  exit 1
fi

# Create or overwrite the lastupdated.json file with the current timestamp
echo -e "${BLUE}Creating or overwriting lastupdated.json with the current timestamp...${NC}"
current_time=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
echo "{\"last_updated\": \"$current_time\"}" > "$LAST_UPDATED_FILE"

# Upload the directories to S3
echo -e "${BLUE}Uploading $LOCAL_DIR to $S3_DIR...${NC}"
if aws s3 cp "$LOCAL_DIR" "$S3_DIR" --recursive; then
  echo -e "${GREEN}Directory upload successful.${NC}"

  # Upload the lastupdated.json file to S3
  echo -e "${BLUE}Uploading lastupdated.json to $S3_DIR...${NC}"
  if aws s3 cp "$LAST_UPDATED_FILE" "$S3_DIR/lastupdated.json"; then
    echo -e "${GREEN}lastupdated.json upload successful.${NC}"
  else
    echo -e "${RED}Failed to upload lastupdated.json.${NC}"
    exit 1
  fi

  # Clean up files in each directory
  echo -e "${BLUE}Cleaning up files in $LOCAL_DIR...${NC}"
  find "$LOCAL_DIR" -type f -exec rm -f {} \;
  echo -e "${GREEN}Files in $LOCAL_DIR have been removed.${NC}"
else
  echo -e "${RED}Upload failed.${NC}"
  exit 1
fi
