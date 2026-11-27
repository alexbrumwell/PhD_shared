#!/bin/bash
#SBATCH --job-name=Het_ANGSD_sharedlist80_pergroup
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

# Combines: (1) Luke's population-level site-discovery step, but now requiring
# data present in >=80% of individuals WITHIN EACH GROUP separately (Modern,
# Historical) and taking the intersection of the two passing site lists --
# mirroring the per-group missingness + bedtools intersect logic used in the
# snpAD/bcftools VCF pipeline -- rather than a single global -minInd cutoff
# applied across the combined cohort; and (2) script 2's per-individual SAF
# step with individual-tailored max-depth cutoffs, looped serially over all
# individuals against the shared site list.

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0 bedtools/2.30.0

# === Inputs ====
LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_samples_selected_6x.list"
AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
FILTERED_META="/projects/erode/people/qhk515/01_OBP/data/OBP_2x_meandepth.list"

MODERN_PREFIX="O-"
HIST_PREFIX="OBP"

MISS_FRAC=0.8   # require data in >=80% of individuals, per group

OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/ANGSD_sharedlist80_pergroup"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

# =========================
# 0. Split combined BAM list into per-group BAM lists
# =========================
# Matches sample names against the same prefixes used in the snpAD pipeline
# (MODERN_PREFIX="O-", HIST_PREFIX="OBP"). Assumes BAM basenames (minus
# .rmdup.bam) start with one of these prefixes, as in script 2's NAME= line.

MODERN_LIST="${OUTDIR}/Modern_bam.list"
HIST_LIST="${OUTDIR}/Historical_bam.list"

# Skip BAM list splitting if already done
if [[ -s "${MODERN_LIST}" && -s "${HIST_LIST}" ]]; then
    echo "$(timestamp) [*] BAM lists already exist, skipping split step"
else
    echo "$(timestamp) [*] Creating per-group BAM lists"
    > "${MODERN_LIST}"
    > "${HIST_LIST}"

    while read -r BAM; do
        NAME=$(basename "${BAM}" .rmdup.bam)
        if [[ "${NAME}" == ${HIST_PREFIX}* ]]; then
            echo "${BAM}" >> "${HIST_LIST}"
        elif [[ "${NAME}" == ${MODERN_PREFIX}* ]]; then
            echo "${BAM}" >> "${MODERN_LIST}"
        else
            echo "WARNING: ${NAME} did not match ${HIST_PREFIX}* or ${MODERN_PREFIX}* -- not assigned to a group" >&2
        fi
    done < "${LIST_COMBINED}"

    n_all=$(wc -l < "${LIST_COMBINED}")
    n_grouped=$(( $(wc -l < "${HIST_LIST}") + $(wc -l < "${MODERN_LIST}") ))
    if [[ "${n_grouped}" -ne "${n_all}" ]]; then
        echo "WARNING: ${n_all} total BAMs but only ${n_grouped} matched a group prefix -- check sample naming" >&2
    fi
fi

[[ -s "${MODERN_LIST}" ]] || { echo "ERROR: ${MODERN_LIST} empty -- no Modern samples found" >&2; exit 1; }
[[ -s "${HIST_LIST}" ]]   || { echo "ERROR: ${HIST_LIST} empty -- no Historical samples found" >&2; exit 1; }

# =========================
# 1. Per-group site discovery, each requiring data in >=80% of THAT group's individuals
# =========================
declare -A GROUP_LIST=( [Modern]="${MODERN_LIST}" [Historical]="${HIST_LIST}" )
declare -A GROUP_SITELIST

for GRP in Modern Historical; do
    LIST="${GROUP_LIST[${GRP}]}"
    N_IND=$(wc -l < "${LIST}")
    MININD=$(awk -v n="${N_IND}" -v p="${MISS_FRAC}" 'BEGIN{v=n*p; iv=int(v); print (v>iv)?iv+1:iv}')
    echo "$(timestamp) [*] ${GRP}: ${N_IND} individuals -> requiring data in >=${MININD} individuals (${MISS_FRAC})"

    # Skip ANGSD if site list already exists
    if [[ -s "${GRP}_80pct.sites.bed" ]]; then
        echo "$(timestamp) [*] ${GRP}_80pct.sites.bed already exists, skipping ANGSD"
    else
        echo "$(timestamp) [*] Starting ANGSD site discovery for ${GRP}"
        angsd -b "${LIST}" -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -out "${GRP}_80pct" \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
            -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd ${MININD} -setMinDepthInd 4 -setMaxDepthInd 50 \
            -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

        echo "$(timestamp) [OK] ANGSD completed for ${GRP}, extracting site list"

        zcat "${GRP}_80pct.mafs.gz" | awk 'NR>1 {print $1"\t"$2-1"\t"$2}' > "${GRP}_80pct.sites.bed"
    fi

    n_sites=$(wc -l < "${GRP}_80pct.sites.bed")
    echo "$(timestamp) [OK] ${GRP}: ${n_sites} sites pass >=${MISS_FRAC} presence"

    GROUP_SITELIST[${GRP}]="${OUTDIR}/${GRP}_80pct.sites.bed"
done

# =========================
# 2. Intersect: sites passing in BOTH groups (this IS the final shared site list)
# =========================
echo "$(timestamp) [*] Intersecting Modern and Historical passing site lists"

# Create a genome file from your reference or scaffolds list (skip if exists)
if [[ ! -f genome.txt ]]; then
    if [[ ! -f ${REF}.fai ]]; then
        samtools faidx ${REF}
    fi
    cut -f1,2 ${REF}.fai > genome.txt
fi

bedtools intersect -sorted -g genome.txt \
    -a "${GROUP_SITELIST[Modern]}" \
    -b "${GROUP_SITELIST[Historical]}" \
    > entire_genome_combined_80pct.bed

awk '{print $1"\t"$3}' entire_genome_combined_80pct.bed > entire_genome_combined_80pct.list

n_shared=$(wc -l < entire_genome_combined_80pct.list)
echo "$(timestamp) [OK] ${n_shared} shared sites pass in both groups"

if [[ "${n_shared}" -eq 0 ]]; then
    echo "WARNING: no shared passing sites across groups" >&2
    exit 1
fi

angsd sites index entire_genome_combined_80pct.list

echo "$(timestamp) [OK] Shared per-group 80%-presence site list ready: ${OUTDIR}/entire_genome_combined_80pct.list"

# === Part 2: per-individual SAF restricted to the shared list, looped serially over ALL individuals (both groups) ===
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

for file in *.ml; do
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
