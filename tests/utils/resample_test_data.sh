#!/bin/sh
# Resample and extract segments from DINO data SSH grid file for use as test data

# Resample to yearly frequency and align on year mid-point (1 July)
cdo yearmean DINO_1m_grid_T.nc DINO_1y_grid_ssh.nc

# Dataset A - extract first two years
cdo seltimestep,1/2 DINO_1y_grid_ssh.nc DINO_subsampled_A/grid_T_2D.nc

# Dataset B (reference dataset) - take years 3 and 4, then shift the year number back to 0001 to match the main set
cdo shifttime,-2year -seltimestep,3/4 DINO_1y_grid_ssh.nc DINO_subsampled_B/grid_T_2D.nc
