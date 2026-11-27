#!/bin/bash 
#SBATCH --job-name=6_concat
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --output=jobs/%x_%j.log
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.22 libdeflate/1.21 htslib vcftools

# === Directories and files ===
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/5_bisnps/miss80"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/6_concat"
mkdir -p $OUTDIR

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

OUTVCF="${OUTDIR}/OBP_merged.26.concat.groupmiss80.bisnps.vcf.gz"


echo "[$(timestamp)] Concatenating scaffold VCFs..."

###############################################################################
# Step 1: Build ordered list of per-scaffold VCFs
###############################################################################

VCF_LIST=()
MISSING=0

while read -r CHRM; do
    # Skip empty lines and comments
    [[ -z "${CHRM}" || "${CHRM}" == \#* ]] && continue

    INVCF="${INDIR}/OBP_merged.26.${CHRM}.groupmiss80.bisnps.vcf.gz"

    if [[ -f "${INVCF}" ]]; then
        VCF_LIST+=("${INVCF}")
    else
        echo "ERROR: Missing VCF: ${INVCF}" >&2
        MISSING=$((MISSING + 1))
    fi
done < "${CHROM_LIST}"


if (( ${#VCF_LIST[@]} == 0 )); then
    echo "ERROR: No input VCFs were found." >&2
    exit 1
fi

echo "[$(timestamp)] Found ${#VCF_LIST[@]} scaffold VCFs."

###############################################################################
# Step 2: Concatenate scaffold VCFs
###############################################################################

bcftools concat \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "${VCF_LIST[@]}" \
    -Oz \
    -o "${OUTVCF}"

###############################################################################
# Step 3: Index and calculate summary statistics
###############################################################################

bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    "${OUTVCF}"

bcftools stats \
    -s - \
    "${OUTVCF}" \
    > "${OUTVCF%.vcf.gz}.stats"

echo "[$(timestamp)] Finished concatenating biallelic SNP VCFs."
echo "Output: ${OUTVCF}"
