#!/usr/bin/env bash
# Build the test-data fixture archive consumed by CI.
#
# Pipeline:
#   1. cdo yearmean monthly 2D file -> annual 2D file
#   2. cdo seltimestep / shifttime to build A (years 1-2) and B (years 3-4
#      relabelled as 1-2) for every time-evolving grid file
#   3. python create_test_data.py to write restart.nc and mesh_mask.nc
#      (with deterministic noise on restart B)
#   4. run spinup_evaluation three times to populate results/
#   5. tar everything into spinup-eval-fixtures.tar.gz
#
# Inputs (paths relative to repo root):
#   data/DINO/DINO_1m_grid_T.nc          monthly 2D source
#   data/DINO/DINO_1y_grid_T.nc          annual 3D T (toce, soce, ...)
#   data/DINO/DINO_1y_grid_U.nc          annual 3D U
#   data/DINO/DINO_1y_grid_V.nc          annual 3D V
#   data/DINO/DINO_00576000_restart.nc   restart
#   data/DINO/mesh_mask.nc               mesh mask
#
# Outputs (under tests/utils/build_fixtures/out/):
#   DINO_subsampled_A/{grid_T_2D,grid_T_3D,grid_U_3D,grid_V_3D,restart,mesh_mask}.nc
#   DINO_subsampled_B/<same>
#   results/{diff,original,self}/
#   spinup-eval-fixtures.tar.gz
#
# Requirements: python (project installed), cdo, tar.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

SRC="data/DINO"
OUT_DIR="tests/utils/build_fixtures_out"
SUB_A="$OUT_DIR/DINO_subsampled_A"
SUB_B="$OUT_DIR/DINO_subsampled_B"
RESULTS="$OUT_DIR/results"
ARCHIVE="$OUT_DIR/spinup-eval-fixtures.tar.gz"

# canonical fixture name -> source file
declare -a CANON=(grid_T_2D.nc grid_T_3D.nc grid_U_3D.nc grid_V_3D.nc)

# --- preflight ---------------------------------------------------------------
for f in DINO_1m_grid_T.nc DINO_1y_grid_T.nc DINO_1y_grid_U.nc DINO_1y_grid_V.nc \
         DINO_00576000_restart.nc mesh_mask.nc; do
    [[ -f "$SRC/$f" ]] || { echo "missing $SRC/$f" >&2; exit 1; }
done
command -v cdo >/dev/null    || { echo "cdo not on PATH" >&2; exit 1; }
command -v python >/dev/null || { echo "python not on PATH" >&2; exit 1; }

[[ -d "$RESULTS" ]] && rm -r "$RESULTS"
[[ -e "$ARCHIVE" ]] && rm "$ARCHIVE"
[[ -d "$SUB_A" ]] && find "$SUB_A" -maxdepth 1 -name '*.nc' -delete
[[ -d "$SUB_B" ]] && find "$SUB_B" -maxdepth 1 -name '*.nc' -delete
mkdir -p "$SUB_A" "$SUB_B" "$RESULTS"

STAGE="$(mktemp -d)"
trap 'rm -r "$STAGE"' EXIT

# --- 1. monthly -> annual for the 2D file ------------------------------------
echo "[1/5] cdo yearmean monthly 2D -> annual"
ANNUAL_2D="$STAGE/grid_T_2D_yearly.nc"
# Force the time axis to mid-year on the 360-day calendar so the 2D file
# matches the NEMO Jul-1 convention used by the 3D files. cftime rejects
# "years" units on 360_day, so use a 360-day increment.
cdo yearmean \
    "$SRC/DINO_1m_grid_T.nc" "$ANNUAL_2D"

# Stage canonical annual sources together
cp "$ANNUAL_2D"             "$STAGE/grid_T_2D.nc"
cp "$SRC/DINO_1y_grid_T.nc" "$STAGE/grid_T_3D.nc"
cp "$SRC/DINO_1y_grid_U.nc" "$STAGE/grid_U_3D.nc"
cp "$SRC/DINO_1y_grid_V.nc" "$STAGE/grid_V_3D.nc"

# --- 2. cdo time slicing for A and B -----------------------------------------
echo "[2/5] cdo seltimestep / shifttime -> A and B"
for fname in "${CANON[@]}"; do
    in="$STAGE/$fname"
    cdo -s seltimestep,1/2 "$in" "$SUB_A/$fname"
    cdo -s shifttime,-2year -seltimestep,3/4 "$in" "$SUB_B/$fname"
done

# --- 3. restart + mesh_mask via python ---------------------------------------
echo "[3/5] python create_test_data.py (restart + mesh_mask)"
# create_test_data.py expects restart.nc and mesh_mask.nc in src_dir, so
# stage symlinks with the canonical names.
mkdir -p "$STAGE/canon_src"
ln -sf "$REPO_ROOT/$SRC/DINO_00576000_restart.nc" "$STAGE/canon_src/restart.nc"
ln -sf "$REPO_ROOT/$SRC/mesh_mask.nc"             "$STAGE/canon_src/mesh_mask.nc"
python tests/utils/create_test_data.py "$STAGE/canon_src" "$SUB_A" "$SUB_B"

# --- 4. run spinup_evaluation ------------------------------------------------
echo "[4/5] running spinup_evaluation (diff/original/self)"
python -m nemo_spinup_evaluation --sim-path="./$SUB_A/" --ref-sim-path="./$SUB_B/" \
    --mode both --results-dir="./$RESULTS/diff"
python -m nemo_spinup_evaluation --sim-path="./$SUB_A/" \
    --mode both --results-dir="./$RESULTS/original"
python -m nemo_spinup_evaluation --sim-path="./$SUB_A/" --ref-sim-path="./$SUB_A/" \
    --mode both --results-dir="./$RESULTS/self"

# --- 5. archive --------------------------------------------------------------
echo "[5/5] creating $ARCHIVE"
tar czf "$ARCHIVE" -C "$OUT_DIR" \
    DINO_subsampled_A DINO_subsampled_B

echo "done: $ARCHIVE"
