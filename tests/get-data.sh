#! /bin/bash
curl -O https://thredds-su.ipsl.fr/thredds/fileServer/NEMO-SPINUP/evaluation/test-data/dino-test.tar.gz
# extract to tests/data
echo "Extracting test data to tests/data"
mkdir -p tests/data
tar -xzf dino-test.tar.gz -C tests/data
