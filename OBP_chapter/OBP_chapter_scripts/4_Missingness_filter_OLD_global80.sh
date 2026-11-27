#!/bin/bash
#SBATCH --job-name=4_Missingness_filter
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --array=1-357%100
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
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/3_MergeGVCFs_perChrm"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_stricter_tests/"
mkdir -p "${OUTDIR}"

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHROM_LIST}")

INVCF="${INDIR}/OBP_merged.26.${chrom}.filt.g.vcf.gz"

if [[ ! -s "${INVCF}" ]]; then
    echo "ERROR: Input VCF not found: ${INVCF}"
    exit 1
fi

echo "$(timestamp) Filtering ${chrom}"

# === 80% called filter ===
# For 26 samples, this means that there should be at least 20 samples 

OUT80="${OUTDIR}/OBP_merged.26.${chrom}.miss80.vcf.gz"

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -i 'F_MISSING<=0.20' \
    -Oz \
    -o "${OUT80}" \
    "${INVCF}"

bcftools index -f "${OUT80}"
bcftools stats -s - "${OUT80}" > "${OUT80%.gz}.stats"


# === 100% called filter ===
# F_MISSING == 0 means no missing genotypes across all samples.

OUT100="${OUTDIR}/OBP_merged.26.${chrom}.miss100.vcf.gz"

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -i 'F_MISSING==0' \
    -Oz \
    -o "${OUT100}" \
    "${INVCF}"

bcftools index -f "${OUT100}"
bcftools stats -s - "${OUT100}" > "${OUT100%.gz}.stats"


# === Quick site counts ===
echo -e "Chromosome\tInput_sites\tMissing90_sites\tMissing100_sites" > "${OUTDIR}/missingness_site_counts.tsv"

echo -e "${chrom}\t$(bcftools view -H "${INVCF}" | wc -l)\t$(bcftools view -H "${OUT80}" | wc -l)\t$(bcftools view -H "${OUT100}" | wc -l)" \
    >> "${OUTDIR}/missingness_site_counts.tsv"

echo "$(timestamp) Done filtering ${chrom}"
echo "$(timestamp) 80% VCF:  ${OUT80}"
echo "$(timestamp) 100% VCF: ${OUT100}"
