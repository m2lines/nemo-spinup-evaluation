#!/bin/bash
curl -O https://thredds-su.ipsl.fr/thredds/fileServer/NEMO-SPINUP/evaluation/test-data/subsampled/current/DINO-simple-small.tar.gz
echo "Extracting test data to tests/data"
tar -xzf DINO-simple-small.tar.gz -C tests/data
rm DINO-simple-small.tar.gz
