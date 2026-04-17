#!/bin/bash
set -e

DATA_DIR="tests/data"
ZIP_FILE="evaluation-test.zip"
ZENODO_URL="https://zenodo.org/records/19557419/files/${ZIP_FILE}"

mkdir -p "${DATA_DIR}"
curl -L -o "${ZIP_FILE}" "${ZENODO_URL}"
unzip -o "${ZIP_FILE}" -d "${DATA_DIR}"
rm "${ZIP_FILE}"
