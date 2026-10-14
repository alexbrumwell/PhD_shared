#!/bin/bash
#SBATCH --job-name=tests_scaf1_MAF_down5x_part1
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0 bcftools/1.22

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/down5x"

mkdir -p ${DIR}

cd ${DIR}

# Variables
# Downsampled list of BAMs
BAM_LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/5x/Down_5x.list"

# Variables
CALLABLE_SITES_LUKE_COMBINED="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het/entire_genome_combined_nosex.list" #total sites 717644389 
LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"

AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_scaffolds_autosomal.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
REF_INDEX="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"

CALLABLE_SITES_MAF_005="maf_filter_sites/scaf1_callable_sites_minmaf0.05"
CALLABLE_SITES_MAF_002="maf_filter_sites/scaf1_callable_sites_minmaf0.02"

# PRIOR TO SUBMITTING JOB ARRAY - subset Lukes callable sites to scaffold 1
#awk '$1=="scaffold_1"{print; next} $1!="scaffold_1"{exit}' ${CALLABLE_SITES_LUKE_COMBINED} > scaf1_callable_sites.txt
#awk 'BEGIN{OFS="\t"} {$1=$1; print}' scaf1_callable_sites.txt > scaf1_callable_sites.txt.tab && angsd sites scaf1_callable_sites.txt.tab

# CREATE NEW LIST OF CALLABLE SITES ON SCAF 1 USING MAF FILTER
mkdir -p maf_filter_sites

# skip if maf already exists and is not empty
if [ -s "${CALLABLE_SITES_MAF_005}.mafs.gz" ]; then
    echo "$(timestamp) [INFO] MAF file for minmaf 0.05 already exists and is not empty. Skipping ANGSD step."
else
  angsd -b "${BAM_LIST_COMBINED}" -ref "${REF}" -anc "${REF}" -r scaffold_1 \
    -sites "${LUKE_LIST_SCAF1}" \
    -out ${CALLABLE_SITES_MAF_005} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd 2 \
    -dumpCounts 3 \
    -doMajorMinor 1 -GL 2 -doGlf 2  -doMaf 2 -minMaf 0.05 \
    -nThreads "${SLURM_CPUS_PER_TASK}"
fi
  zcat ${CALLABLE_SITES_MAF_005}.mafs.gz | awk 'NR>1 {print $1 "\t" $2}' > ${CALLABLE_SITES_MAF_005}.txt

  angsd sites index ${CALLABLE_SITES_MAF_005}.txt

echo "$(timestamp) [OK] Callable sites with minmaf 0.05 created and indexed."

# Now create one with minmaf 0.02
if [ -s "${CALLABLE_SITES_MAF_002}.mafs.gz" ]; then
    echo "$(timestamp) [INFO] MAF file for minmaf 0.02 already exists and is not empty. Skipping ANGSD step."
else
  angsd -b "${BAM_LIST_COMBINED}" -ref "${REF}" -anc "${REF}" -r scaffold_1 \
    -sites "${LUKE_LIST_SCAF1}" \
    -out ${CALLABLE_SITES_MAF_002} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd 2 \
    -dumpCounts 3 \
    -doMajorMinor 1 -GL 2 -doGlf 2  -doMaf 2 -minMaf 0.02 \
    -nThreads "${SLURM_CPUS_PER_TASK}"
fi
  zcat ${CALLABLE_SITES_MAF_002}.mafs.gz | awk 'NR>1 {print $1 "\t" $2}' > ${CALLABLE_SITES_MAF_002}.txt
  
  angsd sites index ${CALLABLE_SITES_MAF_002}.txt

echo "$(timestamp) [OK] Callable sites with minmaf 0.02 created and indexed."

exit 0



#!/bin/bash
#SBATCH --job-name=tests_scaf1_MAF_down5x_part2_0.05
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%10
#SBATCH --output=jobs/%x_%A_%a.out
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk


set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0 bcftools/1.22

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/down5x"

mkdir -p ${DIR}

cd ${DIR}

# Variables
# Downsampled list of BAMs
BAM_LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/5x/Down_5x.list"

# Variables
CALLABLE_SITES_LUKE_COMBINED="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het/entire_genome_combined_nosex.list" #total sites 717644389 
LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"

AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_scaffolds_autosomal.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
REF_INDEX="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"

CALLABLE_SITES_MAF_005="maf_filter_sites/scaf1_callable_sites_minmaf0.05"
CALLABLE_SITES_MAF_002="maf_filter_sites/scaf1_callable_sites_minmaf0.02"


# Run ANGSD per indivudal
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BAM_LIST_COMBINED})
if [ -z "${SAMPLE}" ]; then
    echo "$(timestamp) [ERROR] No sample found for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

NAME=$(basename "${SAMPLE}" .bam)
NAME=${NAME%.down5x}
NAME=${NAME%.rmdup}

echo "$(timestamp) [*] Starting ANGSD for individual ${NAME}"


# Filtered ANGSD - # Requierements (at least >4x coverage and <50x coverage per site to be included in the analysis)
angsd -i ${SAMPLE} -ref ${REF} -r scaffold_1 -anc ${REF} \
    -sites ${CALLABLE_SITES_MAF_005}.txt  -out maf_filter_sites/${NAME}_nosex_filtlist_0.05 \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
    -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 50 -doCounts 1 -dumpCounts 3 \
    -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1 
    

echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

# Now convert calculate SFS and convert to wimSFS format 
winsfs maf_filter_sites/${NAME}_nosex_filtlist_0.05.saf.idx > maf_filter_sites/${NAME}_filter_0.05_wim.ml

echo "$(timestamp) [OK] winSFS completed for ${NAME}"


exit 0

#!/bin/bash
#SBATCH --job-name=tests_scaf1_MAF_part2_0.02
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%10
#SBATCH --output=jobs/%x_%A_%a.out
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk


set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0 bcftools/1.22

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/down5x"

mkdir -p ${DIR}

cd ${DIR}

# Variables
# Downsampled list of BAMs
BAM_LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/5x/Down_5x.list"

# Variables
CALLABLE_SITES_LUKE_COMBINED="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het/entire_genome_combined_nosex.list" #total sites 717644389 
LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"

AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_scaffolds_autosomal.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
REF_INDEX="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"

CALLABLE_SITES_MAF_005="maf_filter_sites/scaf1_callable_sites_minmaf0.05"
CALLABLE_SITES_MAF_002="maf_filter_sites/scaf1_callable_sites_minmaf0.02"

# Run ANGSD per indivudal
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BAM_LIST_COMBINED})
if [ -z "${SAMPLE}" ]; then
    echo "$(timestamp) [ERROR] No sample found for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

NAME=$(basename "${SAMPLE}" .bam)
NAME=${NAME%.down5x}
NAME=${NAME%.rmdup}

echo "$(timestamp) [*] Starting ANGSD for individual ${NAME}"


# Filtered ANGSD - # Requierements (at least >4x coverage and <50x coverage per site to be included in the analysis)
angsd -i ${SAMPLE} -ref ${REF} -r scaffold_1 -anc ${REF} \
    -sites ${CALLABLE_SITES_MAF_002}.txt  -out maf_filter_sites/${NAME}_nosex_filtlist_0.02 \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
    -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 50 -doCounts 1 -dumpCounts 3 \
    -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1 

echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

# Now convert calculate SFS and convert to wimSFS format 
winsfs maf_filter_sites/${NAME}_nosex_filtlist_0.02.saf.idx > maf_filter_sites/${NAME}_filter_0.02_wim.ml

echo "$(timestamp) [OK] winSFS completed for ${NAME}"


exit 0

# AFTER RUNNING SCRIPT

# NOw that I've run the scripts I am interested in counting the hets that have been kept that are higher confidence.
# To do this I will use code to collect the ouput

SUMMARY="heterozygosity_summary_0.05.tsv"

echo -e "sample\ttotal_sites\thet_sites\thet" > "${SUMMARY}"

shopt -s nullglob

for file in *0.05*.ml; do
  sample=$(basename "${file}" _filter_wim.ml)

  grep -v "SHAPE" "${file}" | awk -v sample="${sample}" '
    {
      total_sites += ($1 + $2 + $3)
      het_sites   += $2
    }
    END {
      if (total_sites > 0) {
        het = het_sites / total_sites
        printf("%s\t%.0f\t%.0f\t%.8f\n", sample, total_sites, het_sites, het)
      } else {
        printf("%s\t0\t0\tNA\n", sample)
      }
    }' >> "${SUMMARY}"
done

echo "Wrote: ${SUMMARY}"

# And now that I have the values I can add those to the R dataframe that I created previously 