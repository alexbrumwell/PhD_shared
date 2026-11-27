#!/bin/bash
#SBATCH --job-name=Het_ANGSD_sharedlist80_global
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# Combines: (1) Luke's population-level site-discovery step, but requiring
# data present in >=80% of individuals per site (rather than -minInd 2), and
# (2) script 2's per-individual SAF step with individual-tailored max-depth
# cutoffs, looped serially over all individuals against the shared site list.

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0

# === Inputs ====
LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_35.list"
AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
FILTERED_META="/projects/erode/people/qhk515/01_OBP/data/OBP_2x_meandepth.list"

OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/ANGSD_sharedlist80_global/"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

# === Part 1: build shared site list requiring data in >=80% of individuals ===
N_IND=$(wc -l < "${LIST_COMBINED}")
MININD=$(awk -v n="${N_IND}" -v p=0.8 'BEGIN{v=n*p; iv=int(v); print (v>iv)?iv+1:iv}')
echo "$(timestamp) [*] ${N_IND} individuals in cohort -> requiring data in >=${MININD} individuals per site (80%)"

echo "$(timestamp) [*] Starting ANGSD to build shared 80%-presence site list"
angsd -b ${LIST_COMBINED} -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -out entire_genome_combined_80pct \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd ${MININD} -setMinDepthInd 4 -setMaxDepthInd 50 \
    -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

echo "$(timestamp) [OK] ANGSD completed, extracting site list"

zcat entire_genome_combined_80pct.mafs.gz | awk '{print $1,$2}' | sed 1d > entire_genome_combined_80pct.list && angsd sites index entire_genome_combined_80pct.list



echo "$(timestamp) [OK] Shared 80%-presence site list ready: ${OUTDIR}/entire_genome_combined_80pct.list"

# === Part 2: per-individual SAF restricted to the shared list, looped serially ===
SITES_80PCT="${OUTDIR}/entire_genome_combined_80pct.list"

for BAM in $(cat ${LIST_COMBINED}); do
    NAME=$(basename "${BAM}" .rmdup.bam)
    maxdepth=$(awk -v s="$NAME" '$1 == s {print $3}' "$FILTERED_META")

    echo "$(timestamp) [*] Starting ANGSD for individual ${NAME} (maxdepth=${maxdepth})"

    angsd -i ${BAM} -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -sites ${SITES_80PCT} -out ${NAME}_sharedlist80_min4x \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
        -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth ${maxdepth} -doCounts 1 \
        -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1

    echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

    winsfs ${NAME}_sharedlist80_min4x.saf.idx > ${NAME}_sharedlist80_min4x.ml

    echo "$(timestamp) [OK] winSFS completed for ${NAME}"
done



# To get summary het values across all individuals, run the following code after the job is done:

SUMMARY="heterozygosity_summary.tsv"

echo -e "sample\tn_boot\tmean_H\tsd_H\tmean_sites\tsd_sites" > "${SUMMARY}"

shopt -s nullglob

for file in *_sharedlist80_min4x.ml; do
  sample=$(basename "${file}" _sharedlist80_min4x.ml)

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

