#!/bin/bash
#SBATCH --job-name=3_Merge
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --array=1
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21
module load libdeflate/1.21 htslib

timestamp() { date '+%F %T'; }

# === Directories ===

INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/2_FiltGVCFs"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/3_MergeGVCFs_perChrm"
mkdir -p ${OUTDIR}

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CHROM_LIST )


echo "$(timestamp) Merging gVCFs for $chrom ..."

# === Create VCF list ===
VCF_LIST="${OUTDIR}/vcfs_merge.${chrom}.list"
: > "${VCF_LIST}"

missing=0

while IFS= read -r sample || [[ -n "${sample}" ]]; do
    [[ -z "${sample}" ]] && continue

    v="${INDIR}/${chrom}/${sample}_${chrom}.filt.gvcf.gz"

    if [[ ! -s "${v}" ]]; then
        echo "WARNING: Missing VCF: ${v}"
        missing=$((missing + 1))
        continue
    fi

    if [[ ! -s "${v}.csi" && ! -s "${v}.tbi" ]]; then
        echo "$(timestamp) Indexing missing index for ${v}"
        bcftools index -f "${v}"
    fi

    echo "${v}" >> "${VCF_LIST}"

done < "${LIST_SELECTED_SAMPLES}"

n_vcfs=$(wc -l < "${VCF_LIST}")

if [[ "${n_vcfs}" -eq 0 ]]; then
    echo "ERROR: No valid VCFs found for ${chrom}"
    exit 1
fi

echo "$(timestamp) Found ${n_vcfs} VCFs for ${chrom}; missing ${missing}"

# === Merge ===
OUTVCF="${OUTDIR}/OBP_merged.26.${chrom}.filt.g.vcf.gz"

bcftools merge \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -l "${VCF_LIST}" \
    -Ou \
| bcftools +fill-tags \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "${OUTVCF}" \
    -- -t AC,AN,MAF,F_MISSING

bcftools index -f "${OUTVCF}"
bcftools stats -s - "${OUTVCF}" > "${OUTVCF%.gz}.stats"

echo "$(timestamp) Done: ${OUTVCF}"


