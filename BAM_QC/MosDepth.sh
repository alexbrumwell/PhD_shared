#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --job-name=Mosdepth_cov
#SBATCH --time=24:00:00

# === SET TARGET DIRECTORY ===
DIR="/projects/erode/people/qhk515/01_OBP/data/bams_outgroups/" 
;lknEKOQN


# === Load modules ===
module load mosdepth/0.3.10   

# === Create output directory ===
OUTDIR="${DIR}/mosdepth_out"
mkdir -p "${OUTDIR}"

# Summary table
SUMMARY_TSV="${OUTDIR}/OBP_mosdepth_mean_coverage.tsv"
echo -e "Sample\tMeanCoverage" > "${SUMMARY_TSV}"

# Optional thresholds (to also get breadth at ≥1×, ≥5×, ≥10×)
THRESH="1,5,10"

for BAM in "${DIR}"*.bam; do
    
    SAMPLE=$(basename "${BAM}" .rmdup.bam)
    OUTFILE="${OUTDIR}/${SAMPLE}"
    
    echo "Running mosdepth for ${SAMPLE}"
    mosdepth --no-per-base --by 100000 --thresholds 1,5,10 "${OUTFILE}" "${BAM}"


    # Parse mean coverage from the summary (the 'total' row)
    MEAN=$(awk '$1=="total"{print $4}' "${OUTFILE}.mosdepth.summary.txt")
    echo -e "${SAMPLE}\t${MEAN}" >> "${SUMMARY_TSV}"
done

echo "Done. Summary table: ${SUMMARY_TSV}"
