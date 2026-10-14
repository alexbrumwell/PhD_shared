#!/bin/bash
#SBATCH --job-name=tests_scaf1_down5x
#SBATCH --time=02:00:00
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

# Downsampled list of BAMs
BAM_LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/5x/Down_5x.list"
INDIR_GVCFS="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down5x/1_VCFs"

# Variables
CALLABLE_SITES_LUKE_COMBINED="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het/entire_genome_combined_nosex.list" #total sites 717644389 
LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"

AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_scaffolds_autosomal.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
REF_INDEX="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"

# PRIOR TO SUBMITTING JOB ARRAY - subset Lukes callable sites to scaffold 1
#awk '$1=="scaffold_1"{print; next} $1!="scaffold_1"{exit}' ${CALLABLE_SITES_LUKE_COMBINED} > scaf1_callable_sites.txt
#awk 'BEGIN{OFS="\t"} {$1=$1; print}' scaf1_callable_sites.txt > scaf1_callable_sites.txt.tab && angsd sites scaf1_callable_sites.txt.tab

# First count the total number of sites in scaffold 1 from the reference genome
RAW_SCAF1_COUNT=$(awk '$1 == "scaffold_1" {print $2}' ${REF_INDEX})

# Count number of callable sites in Lukes list
LUKE_SCAF1_COUNT=$(wc -l < "${LUKE_LIST_SCAF1}")


# For each indivudal, I want to calculate the number of total sites and het sites that remain after I run two methods: ANGSD and snpAD

## METHOD 1: ANGSD
mkdir -p output_angsd

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
    -sites ${LUKE_LIST_SCAF1} -out output_angsd/${NAME}_nosex_filtlist \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
    -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 50 -doCounts 1 -dumpCounts 3 \
    -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1 

echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

# Now convert calculate SFS and convert to wimSFS format 
winsfs output_angsd/${NAME}_nosex_filtlist.saf.idx > output_angsd/${NAME}_filter_wim.ml

echo "$(timestamp) [OK] winSFS completed for ${NAME}"

# Count how many sites are included in the analysis for this individual (i.e. how many sites have coverage between 4x and 50x)
#ANGSD_RAW_SCAF1_COUNT=??????  # I will count ANGSD_RAW_SCAF1_COUNT manuall after ANGSD has finished running by going through the log file and looking at the "Total number of sites analyzed"

ANGSD_TOTAL_SITES=$(cat output_angsd/${NAME}_filter_wim.ml | awk 'NR>1  {print $1 + $2 + $3}')
ANGSD_HET_SITES=$(cat output_angsd/${NAME}_filter_wim.ml | awk 'NR>1 {print $2}')

# Previous alternative method - Compare count to ANGSD genotype calls
#ANGSD_TOTAL_SITES=$(awk '{if($1 != -1) n++} END{print n+0}' output_angsd/${NAME}_geno.geno)
#ANGSD_HET_SITES=$(awk '{if($1 == 1) n++} END{print n+0}' output_angsd/${NAME}_geno.geno)



## METHOD 2: snpAD
cd ${DIR}
mkdir -p output_snpad

# Get the VCF 
GVCF="${INDIR_GVCFS}/scaffold_1/${NAME}_scaffold_1.vcf.gz"
OUTVCF="output_snpad/${NAME}_scaffold_1_filtered.gvcf.gz"

# Count how many sites are in the raw gVCF in scaf 1
SNPAD_RAW_SCAF1_COUNT=$(bcftools view -H ${GVCF} | wc -l) 

echo "$(timestamp) Starting to filter ${GVCF} "

# Filter gVCF
bcftools view \
  -T "${LUKE_LIST_SCAF1}" \
  -i 'FMT/DP>=4 && FMT/DP<=50' \
  --threads "${SLURM_CPUS_PER_TASK}" \
  -Oz -o "${OUTVCF}" \
  "${GVCF}"

echo "Finished filtering ${GVCF}"

bcftools index -f "${OUTVCF}"

bcftools stats -s - "${OUTVCF}" > "${OUTVCF%.gvcf.gz}.stats"

echo "$(timestamp) Finished stats for ${GVCF}"

# Count how many sites are included in the analysis for this individual from the stats file
SNPAD_TOTAL_SITES=$(grep "PSC" "${OUTVCF%.gvcf.gz}.stats" | tail -n 1 | cut -f 4,5,6 | awk '{print $1 + $2 + $3}')
SNPAD_HET_SITES=$(grep "PSC" "${OUTVCF%.gvcf.gz}.stats" | tail -n 1 | cut -f 6)

#Alternative count methods (although I've noticed that it tends to overinflate het counts):
#SNPAD_TOTAL_SITES=$(bcftools query -f '[%GT\n]' "${OUTVCF}" | awk '$1 != "./." && $1 != "." {n++} END{print n+0}')
#SNPAD_HET_SITES=$(bcftools query -f '[%GT\n]' "${OUTVCF}" | awk '$1=="0/1" || $1=="1/0" || $1=="0|1" || $1=="1|0" {n++} END{print n+0}')



# Print all results in one line for easy comparison
SUMMARY_SITES_DOC="summary_sites_scaf1.txt"
# print header only once
if [ ! -f "${SUMMARY_SITES_DOC}" ]; then
    echo -e "Individual\tRaw_Scaf1_Count\tLuke_Scaf1_Count\tANGSD_Raw_Scaf1_Count\tANGSD_Total_Sites\tANGSD_Het_Sites\tSNPAD_Raw_Scaf1_Count\tSNPAD_Total_Sites\tSNPAD_Het_Sites" > ${SUMMARY_SITES_DOC}
fi
echo -e "${NAME}\t${RAW_SCAF1_COUNT:-NA}\t${LUKE_SCAF1_COUNT:-NA}\t${ANGSD_RAW_SCAF1_COUNT:-NA}\t${ANGSD_TOTAL_SITES:-NA}\t${ANGSD_HET_SITES:-NA}\t${SNPAD_RAW_SCAF1_COUNT:-NA}\t${SNPAD_TOTAL_SITES:-NA}\t${SNPAD_HET_SITES:-NA}" >> ${SUMMARY_SITES_DOC}

echo "$(timestamp) Finished printing summary values to ${SUMMARY_SITES_DOC}"


exit 0

# See R script after