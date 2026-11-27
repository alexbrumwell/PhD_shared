#!/bin/bash
#SBATCH --job-name=2_Filter
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --array=1-357
#SBATCH --time=02:00:00
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

module load perl/5.38.0 gsl/2.5 bcftools/1.21

timestamp() { date '+%F %T'; }


INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/1_VCFs"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/2_FiltGVCFs"
FILTERED_META="/projects/erode/people/qhk515/01_OBP/data/OBP_2x_meandepth.list"

TORUN=$1

read -r chrom sample < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TORUN")

maxdepth=$(awk -v s="$sample" '$1 == s {print $3}' "$FILTERED_META")

mkdir -p "${OUTDIR}/${chrom}"

INVCF="${INDIR}/${chrom}/${sample}_${chrom}.vcf.gz"
OUTVCF="${OUTDIR}/${chrom}/${sample}_${chrom}.filt.gvcf.gz"
TMPVCF="${OUTDIR}/${chrom}/${sample}_${chrom}.tmp.vcf.gz"
SAMPLEFILE="${OUTDIR}/${chrom}/${sample}_${chrom}.sample.txt"

echo "$(timestamp)"
echo "Filtering ${sample} ${chrom}"
echo "Input: $INVCF"
echo "Output: $OUTVCF"
echo "Max depth: ${maxdepth}"


bcftools +setGT "${INVCF}" -Ou -- -t q -n '.' -i "FMT/DP<4 || FMT/GQ<30 || FMT/DP>${maxdepth} " | bcftools view -Oz -o "$TMPVCF"
#bcftools filter -e "FMT/DP<4 || FMT/GQ<30 || FMT/DP>${maxdepth} " "${INVCF}" -Oz -o "$TMPVCF"

bcftools index -f "$TMPVCF"

# Now reheader file (usually snpad names files correctly but I messed up at an earlier stage so I have to rename them here)
echo "$(timestamp) Reheadering"

echo "$sample" > "$SAMPLEFILE"

bcftools reheader -s $SAMPLEFILE -o "$OUTVCF" "$TMPVCF"
bcftools index -f "$OUTVCF"


bcftools stats -s - "$OUTVCF" > "${OUTVCF%.gvcf.gz}.stats"

rm -f "$TMPVCF" "$TMPVCF.csi" "$SAMPLEFILE"

echo "$(timestamp) Finished ${sample} ${chrom}"
