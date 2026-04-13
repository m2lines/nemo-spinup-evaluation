#!/bin/bash
set -e

DATA_DIR="tests/data"
ZIP_FILE="evaluation-test.zip"
CONCEPT_ID="19474413"

# Resolve the concept ID to the latest version
LATEST=$(curl -sI "https://zenodo.org/records/${CONCEPT_ID}" | grep -i "^location:" | tr -d '\r' | awk -F/ '{print $NF}')
ZENODO_URL="https://zenodo.org/records/${LATEST}/files/${ZIP_FILE}"

echo $LATEST

mkdir -p "${DATA_DIR}"
curl -L -o "${ZIP_FILE}" "${ZENODO_URL}"
unzip -o "${ZIP_FILE}" -d "${DATA_DIR}"
# rm "${ZIP_FILE}"
