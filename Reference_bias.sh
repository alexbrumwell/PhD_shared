#!/bin/bash
#SBATCH --job-name=IBS_refbias
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err

set -euo pipefail

# ====== MODULES ======
module purge
module load angsd/0.940

mkdir -p jobs

# ====== USER INPUTS (EDIT THESE) ======
OUT_DIR="/projects/erode/people/qhk515/01_OBP/hdna_tests/reference_bias/"
OUT_PREFIX="${OUT_DIR}OBP_IBS_filtered" # Prefix for output files
BAM_LIST="/projects/erode/people/qhk515/01_OBP/hdna_tests/reference_bias/scripts/OBPs_all.list"
COMMON_SITES="/projects/erode/people/qhk515/01_OBP/data/angsd_callable/all/historical/miss0/historical_all.sites" #specify regions good quality

# ===== Other variables ======
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
mkdir -p "${OUT_DIR}"

THREADS=8


# ====== TEST 1: IBS without list of sites ===========================
echo "[$(date)] Starting ANGSD IBS run (no MQ/BQ filters) on all sites..."
angsd -bam "${BAM_LIST}" -ref "${REF}" \
      -doCounts 1 \
      -doIBS 1 \
      -makeMatrix 1 \
      -minInd 1 \
      -P ${THREADS} \
      -out "${OUT_PREFIX}"

echo "[$(date)] Done. Files:"
echo "  IBS distance matrix : ${OUT_PREFIX}.ibsMat"
echo "  Sampled bases       : ${OUT_PREFIX}.ibs.gz (single-read per site)"

exit 0


# ====== TEST 2: IBS WITH list of sites ===========================
echo "[$(date)] Starting ANGSD IBS run (no MQ/BQ filters) with list of sites..."
angsd -bam "${BAM_LIST}" -ref "${REF}" \
      -sites "${COMMON_SITES}" \
      -doCounts 1 \
      -doIBS 1 \
      -makeMatrix 1 \
      -minInd 1 \
      -P ${THREADS} \
      -out "${OUT_PREFIX}"

echo "[$(date)] Done. Files:"
echo "  IBS distance matrix : ${OUT_PREFIX}.ibsMat"
echo "  Sampled bases       : ${OUT_PREFIX}.ibs.gz (single-read per site)"

exit 0



# ------------------------------------------------------------------------------
# LATER 

# Prefix used in ANGSD IBS run (without extension)
OUT_DIR="/projects/erode/people/qhk515/01_OBP/hdna_tests/reference_bias/"
OUT_PREFIX="${OUT_DIR}OBP_IBS_filtered" # Prefix for output files

# List of BAMs used in ANGSD, in the same order
BAM_LIST="/projects/erode/people/qhk515/01_OBP/hdna_tests/reference_bias/scripts/OBPs_all.list"

# Output matrix (labelled)
OUT_TSV="${OUT_PREFIX}_matrix.tsv"

IBSMAT="${OUT_PREFIX}.ibsMat"

echo "Processing IBS matrix: ${IBSMAT}"
echo "Using BAM list       : ${BAM_LIST}"
echo "Writing TSV          : ${OUT_TSV}"

# Get sample names from BAM list
awk -F/ '{gsub(/\.(rmdup\.)?bam$/,"",$NF); print $NF}' "${BAM_LIST}" > "${OUT_PREFIX}.samplenames"

# Read sample names into an array
mapfile -t SAMPLES < "${OUT_PREFIX}.samplenames"

# Write header + matrix with row names
{
  # header
  printf "Sample"
  for s in "${SAMPLES[@]}"; do
    printf "\t%s" "${s}"
  done
  printf "\n"

  # rows: prepend row name to each line of the ibsMat
  i=0
  while IFS= read -r line; do
    printf "%s\t%s\n" "${SAMPLES[$i]}" "${line}"
    ((i++))
  done < "${IBSMAT}"

} > "${OUT_TSV}"

echo "Done. Labelled IBS matrix written to: ${OUT_TSV}"