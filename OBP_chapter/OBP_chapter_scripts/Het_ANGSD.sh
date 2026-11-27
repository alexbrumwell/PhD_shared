#!/bin/bash
#SBATCH --job-name=Het_ANGSD
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%5
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0

# === Inputs ====
LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_35.list"
AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
FILTERED_META="/projects/erode/people/qhk515/01_OBP/data/OBP_2x_meandepth.list"

OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/het_ANGSD_ind"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST_COMBINED}")
NAME=$(basename "${BAM}" .rmdup.bam)
maxdepth=$(awk -v s="$NAME" '$1 == s {print $3}' "$FILTERED_META")

# Step 1: Calculate high quality sites per individual 
echo "$(timestamp) [*] Creating sites list for individual ${NAME}"

angsd -i ${BAM} -ref ${REF} -rf ${AUTOSOMAL_SCAFFOLDS} -anc ${REF} \
-out ${NAME}_rmtrans_GL2_min4x \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans 1 -trim 0 -C 50 \
-baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth $maxdepth -doCounts 1 -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 \
-nThreads ${SLURM_CPUS_PER_TASK}

zcat ${NAME}_rmtrans_GL2_min4x.mafs.gz| awk '{print $1,$2}' | sed 1d > ${NAME}_rmtrans_GL2_min4x.list && angsd sites index ${NAME}_rmtrans_GL2_min4x.list

# Step 2: Run ANGSD
echo "$(timestamp) [*] Creating GL for individual ${NAME}"

angsd  -i ${BAM} -ref $REF -sites ${NAME}_rmtrans_GL2_min4x.list -anc $REF -out ${NAME}_rmtrans_GL2_min4x_FINAL -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -noTrans 1 -C 50 -baq 0 -minMapQ 30 -minQ 20  -setMinDepth 4 -setMaxDepth $maxdepth -doCounts 1 -nThreads $SLURM_CPUS_PER_TASK -GL 2 -doSaf 1

# Now convert calculate SFS and convert to wimSFS format 
winsfs ${NAME}_rmtrans_GL2_min4x_FINAL.saf.idx > ${NAME}_rmtrans_GL2_min4x_FINAL.ml

echo "$(timestamp) [OK] winSFS completed for ${NAME}"


exit 0

# ------------------------
# To get summary het values across all individuals, run the following code after all jobs are done:

SUMMARY="heterozygosity_summary.tsv"
:> $SUMMARY

echo -e "sample\tn_boot\tmean_H\tsd_H\tmean_sites\tsd_sites" > "${SUMMARY}"

shopt -s nullglob

for file in *ml; do
  sample=$(basename "${file}" _filter_wim.ml)
  
  grep -v "SHAPE" "${file}" | awk -v sample="${sample}" '
    {
      tot = $1 + $2 + $3
      if (tot > 0) {
        h = $2 / tot
        n++

        sumH   += h
        sumsqH += h*h

        sumT   += tot
        sumsqT += tot*tot
      }
    }
    END {
      if (n == 0) {
        printf("%s\t0\tNA\tNA\tNA\tNA\n", sample)
        exit
      }

      meanH = sumH / n
      if (n > 1) {
        varH = (sumsqH - n*meanH*meanH) / (n-1)
        if (varH < 0) varH = 0
        sdH = sqrt(varH)
      } else {
        sdH = 0
      }

      meanT = sumT / n
      if (n > 1) {
        varT = (sumsqT - n*meanT*meanT) / (n-1)
        if (varT < 0) varT = 0
        sdT = sqrt(varT)
      } else {
        sdT = 0
      }

      printf("%s\t%d\t%.10f\t%.10f\t%.0f\t%.0f\n",
             sample, n, meanH, sdH, meanT, sdT)
    }' >> "${SUMMARY}"
done

echo "Wrote: ${SUMMARY}"
