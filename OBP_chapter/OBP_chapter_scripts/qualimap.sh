# Run Qualimap independently

#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G 
#SBATCH --time=24:00:00
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err

set -euo pipefail

module purge
module load qualimap/2.3 java-jdk/8.0.112
 
 timestamp() { date '+%F %T'; }

# Set directory
DIR=$1

cd $DIR

shopt -s nullglob

OUTDIR="qualimap_reports"
mkdir -p "${OUTDIR}" jobs

for file in ${DIR}/*.bam; do
  sample=$(basename "${file}" .bam)
  
  REPORT_DIR="${OUTDIR}/${sample}_qualimap"

  echo "[$(timestamp)] Running Qualimap for BAM: ${file}"

# Run Qualimap (skip if already done)
  qualimap bamqc \
    -bam "${file}" \
    -outdir "${REPORT_DIR}" \
    -nt "${SLURM_CPUS_PER_TASK}" \
    --java-mem-size=16G

echo "[$(date)] Done for $file"

done 

echo "[$(timestamp)] All Qualimap jobs complete"

exit 0


#Then run MultiQC to aggregate results
module purge
module load qualimap/2.3 java-jdk/8.0.112

 
INDIR="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/qualimap_reports"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/qualimap_reports/multiqc"
mkdir -p $OUTDIR 

# Create a list of the bam qc results (three columns: sample, path, group)
paste <(basename -a *qualimap | sed 's/_qualimap//') <(realpath *qualimap) > ${INDIR}/bam.list
awk 'BEGIN{OFS="\t"} /^O-/  {print $1, $2, "Modern"; next} /^OBP/ {print $1, $2, "Historical"}' ${INDIR}/bam.list > ${INDIR}/bam_group.list

qualimap multi-bamqc \
--data ${INDIR}/bam_group.list \
-outdir $OUTDIR \
-outfile multiqc_OBP \
-outformat PDF:HTML

# Or you can also manually summarise the QC reports

# RUN AFTER:

REPORT_DIR="."
OUT_TSV="${REPORT_DIR}/qualimap_summary.tsv"

echo -e "SampleID\tTotalReads\tMappedReads\tMappingQuality\tMeanCoverage\tStdDevCoverage\tDuplicationRate\tMeanInsertSize\tPctGenome_2x\tPctGenome_5x\tPctGenome_10x" > "${OUT_TSV}"

for report in "${REPORT_DIR}"/*_qualimap/; do
  STATS="${report}/genome_results.txt"
  [[ ! -s "${STATS}" ]] && continue

  SAMPLE=$(basename "${report}" | sed 's/_qualimap$//')

  TOTAL_READS=$(grep "number of reads =" "${STATS}" | awk -F' = ' '{print $2}' | tr -d ',')
  MAPPED_READS=$(grep "number of mapped reads =" "${STATS}" | awk -F' = ' '{print $2}' | tr -d ',')
  MAPQ=$(grep "mean mapping quality =" "${STATS}" | awk -F' = ' '{print $2}')
  MEAN_COV=$(grep "mean coverageData" "${STATS}" | awk -F' = ' '{print $2}')
  STD_COV=$(grep "std coverageData" "${STATS}" | awk -F' = ' '{print $2}')
  DUP_RATE=$(grep "duplication rate" "${STATS}" | awk -F' = ' '{print $2}' | tr -d '% ')
  INSERT_SIZE=$(grep "mean insert size =" "${STATS}" | awk -F' = ' '{print $2}')

  COV2=$(grep "coverageData >= 2X" "${STATS}" | awk '{print $4}' | tr -d '%')
  COV5=$(grep "coverageData >= 5X" "${STATS}" | awk '{print $4}' | tr -d '%')
  COV10=$(grep "coverageData >= 10X" "${STATS}" | awk '{print $4}' | tr -d '%')

  echo -e "${SAMPLE}\t${TOTAL_READS:-NA}\t${MAPPED_READS:-NA}\t${MAPQ:-NA}\t${MEAN_COV:-NA}\t${STD_COV:-NA}\t${DUP_RATE:-NA}\t${INSERT_SIZE:-NA}\t${COV2:-NA}\t${COV5:-NA}\t${COV10:-NA}" \
    >> "${OUT_TSV}"
done

echo "[DONE] Wrote ${OUT_TSV}"


# Later in R -> see R script on local device
