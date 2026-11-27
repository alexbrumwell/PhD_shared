#!/bin/bash
#SBATCH --job-name=4_Missingness_filter
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --array=1-357%100  # set to number of lines in CHROM_LIST
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21 bedtools/2.30.0

timestamp() { date '+%F %T'; }

# =========================
# Parameters
# =========================
MISS_THRESH=0.2   # F_MISSING > this fails, i.e. requires >=80% call rate per group

# =========================
# Directories / inputs
# =========================
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/3_MergeGVCFs_perChrm"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_GroupMissFiltered"
TMPBASE="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_GroupMissFiltered/tmp"
BED_DIR="${OUTDIR}/beds"

CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

MODERN_PREFIX="O-"
HIST_PREFIX="OBP"

LISTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_GroupMissFiltered"
MODERN_LIST="${LISTDIR}/Modern.txt"
HIST_LIST="${LISTDIR}/Historical.txt"

mkdir -p "${OUTDIR}" "${TMPBASE}" "${LISTDIR}" "${BED_DIR}"

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHROM_LIST}")
[[ -n "${chrom}" ]] || { echo "No chrom for array index ${SLURM_ARRAY_TASK_ID}"; exit 0; }

IN_VCF="${INDIR}/OBP_merged.26.${chrom}.filt.g.vcf.gz"
[[ -s "${IN_VCF}" ]] || { echo "ERROR: missing input ${IN_VCF}" >&2; exit 1; }

# =========================
# 0. RUN BEFORE SUBMITTING SCRIPT
# =========================
#        echo "$(timestamp) Generating group lists from ${IN_VCF} sample names"
#        bcftools query -l "${IN_VCF}" > "${LISTDIR}/all_samples.txt"
#        grep "^${HIST_PREFIX}"   "${LISTDIR}/all_samples.txt" > "${HIST_LIST}"
#        grep "^${MODERN_PREFIX}" "${LISTDIR}/all_samples.txt" > "${MODERN_LIST}"

#        n_all=$(wc -l < "${LISTDIR}/all_samples.txt")
#        n_grouped=$(( $(wc -l < "${HIST_LIST}") + $(wc -l < "${MODERN_LIST}") ))
#        if [[ "${n_grouped}" -ne "${n_all}" ]]; then
#            echo "WARNING: ${n_all} total samples but only ${n_grouped} matched a prefix -- check sample naming:" >&2
#            comm -23 <(sort "${LISTDIR}/all_samples.txt") <(cat "${HIST_LIST}" "${MODERN_LIST}" | sort)
#        fi

[[ -s "${MODERN_LIST}" ]] || { echo "ERROR: ${MODERN_LIST} empty/missing after generation" >&2; exit 1; }
[[ -s "${HIST_LIST}" ]]   || { echo "ERROR: ${HIST_LIST} empty/missing after generation" >&2; exit 1; }

TMPDIR="${TMPBASE}/${chrom}"
rm -rf "${TMPDIR}"
mkdir -p "${TMPDIR}"

OUT_VCF="${OUTDIR}/OBP_merged.26.${chrom}.groupmiss80.vcf.gz"
OUT_BED="${BED_DIR}/${chrom}.shared_groupmiss80.bed"

echo "$(timestamp) chrom=${chrom} MISS_THRESH=${MISS_THRESH}"
echo "IN_VCF : ${IN_VCF}"
echo "OUT_VCF: ${OUT_VCF}"



# =========================
# 1. Per-group F_MISSING -> passing sites written straight to a merged BED
# =========================
for GRP in Modern Historical; do
    if [[ "${GRP}" == "Modern" ]]; then LIST="${MODERN_LIST}"; else LIST="${HIST_LIST}"; fi

    bcftools view -S "${LIST}" --force-samples -Ou "${IN_VCF}" \
      | bcftools +fill-tags --threads "${SLURM_CPUS_PER_TASK}" -Ou -- -t F_MISSING \
      | bcftools view -i "F_MISSING<= ${MISS_THRESH}" \
      | bcftools query -f '%CHROM\t%POS0\t%POS\n' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "${TMPDIR}/${GRP}.pass.bed"

    n_int=$(wc -l < "${TMPDIR}/${GRP}.pass.bed")
    echo "$(timestamp) ${GRP}: ${n_int} merged intervals pass <= ${MISS_THRESH} missingness"
done

# =========================
# 2. Intersect: regions passing in BOTH groups (this IS the final mask)
# =========================
bedtools intersect \
    -a "${TMPDIR}/Historical.pass.bed" \
    -b "${TMPDIR}/Modern.pass.bed" \
    > "${OUT_BED}"

n_shared=$(wc -l < "${OUT_BED}")
echo "$(timestamp) ${n_shared} shared intervals pass in both groups -> ${OUT_BED}"

if [[ "${n_shared}" -eq 0 ]]; then
    echo "WARNING: no shared passing regions for ${chrom}"
exit 0
fi

# =========================
# 3. Apply shared BED mask back to the full, all-sample VCF
#    (bcftools recognises the .bed extension as 0-based, half-open)
# =========================
bcftools view -T "${OUT_BED}" -Oz -o "${OUT_VCF}" "${IN_VCF}"
bcftools index -f "${OUT_VCF}"
bcftools stats -s - "${OUT_VCF}" > "${OUT_VCF%.gz}.stats"


echo "$(timestamp) Done: ${OUT_VCF}"
