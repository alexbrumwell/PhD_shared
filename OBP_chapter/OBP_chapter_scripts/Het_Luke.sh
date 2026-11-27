#!/bin/bash
#SBATCH --job-name=Het_Luke
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0

# === Inputs ====
LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_35.list"
LIST_HISTORICAL="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_hist_16.list"
LIST_MODERN="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_mod_19.list"
AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes.list" # All scaffolds except scaffold_5 (sex chromosome)
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"

OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

# Making group lists (requieres at least data from 2 individuals to be included)
# Luke actually only used the combined list for the heterozygosity analysis
echo "$(timestamp) [*] Starting ANGSD for combined groups"
angsd -b ${LIST_COMBINED} -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -out entire_genome_combined_nosex \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd 2 \
    -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK} 

#echo "$(timestamp) [OK] ANGSD completed for combined group, starting modern and historical groups"
#angsd -b ${LIST_MODERN} -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -out entire_genome_modern_nosex \
#    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
#    -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd 2 \
#    -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

#echo "$(timestamp) [OK] ANGSD completed for modern group, starting historical group"
#angsd -b ${LIST_HISTORICAL} -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -out entire_genome_historical_nosex \
#    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
#    -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd 2 \
#   -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

echo "$(timestamp) [OK] Group lists created"

zcat entire_genome_combined_nosex.mafs.gz| awk '{print $1,$2}' | sed 1d > entire_genome_combined_nosex.list && angsd sites index entire_genome_combined_nosex.list
#zcat entire_genome_modern_nosex.mafs.gz| awk '{print $1,$2}' | sed 1d > entire_genome_modern_nosex.list && angsd sites index entire_genome_modern_nosex.list
#zcat entire_genome_historical_nosex.mafs.gz| awk '{print $1,$2}' | sed 1d > entire_genome_historical_nosex.list && angsd sites index entire_genome_historical_nosex.list

echo "$(timestamp) [OK] Group lists filtered to just positions created"


# Run ANGSD per indivudal
for BAM in $(cat ${LIST_COMBINED}); do
    NAME=$(basename "${BAM}" .rmdup.bam)
    echo "$(timestamp) [*] Starting ANGSD for individual ${NAME}"

    # Requierements (at least >4x coverage and <50x coverage per site to be included in the analysis)
    angsd -i ${BAM} -ref ${REF} -rf ${AUTOSOMAL_SCAFFOLDS} -anc ${REF} \
        -sites entire_genome_combined_nosex.list -out ${NAME}_nosex_filtlist \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
        -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 50 -doCounts 1 -dumpCounts 3 \
        -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1

    echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

    # Now convert calculate SFS and convert to wimSFS format 
    winsfs ${NAME}_nosex_filtlist.saf.idx > ${NAME}_filter_wim.ml

    echo "$(timestamp) [OK] winSFS completed for ${NAME}"
done

exit 0
