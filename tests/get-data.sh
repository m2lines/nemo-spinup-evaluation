#! /bin/bash
curl -O https://thredds-su.ipsl.fr/thredds/fileServer/NEMO-SPINUP/evaluation/test-data/subsampled/DINO-simple-small.tar.gz
# extract to tests/data
echo "Extracting test data to tests/data"
mkdir -p tests/data
tar -xzf DINO-simple-small.tar.gz -C tests/data
