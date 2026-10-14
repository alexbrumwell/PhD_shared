#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --job-name=mapDamage_OBPs
#SBATCH --output=/projects/erode/people/qhk515/01_OBP/scripts/OBPs/jobs/mapDamage_OBPs_%j.out
#SBATCH --error=/projects/erode/people/qhk515/01_OBP/scripts/OBPs/jobs/mapDamage_OBPs_%j.err
#SBATCH --time=72:00:00

# Load required modules
module load mapdamage2/2.2.2

# Set variables
LIST="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_OBPs.txt"
DIR_OUTPUT="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/mapdamage_reports/"
REFERENCE="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"  # Reference genome


# Create output directory if it doesn't exist
mkdir -p "$DIR_OUTPUT"

# Loop through each line of the list file
while read -r BAM ID; do
    REPORT_DIR="${DIR_OUTPUT}${ID}_mapdamage"

    # Run mapDamage
    mapDamage -i "$BAM" -r "$REFERENCE" -d "$REPORT_DIR"  --merge-reference-sequences || echo "mapDamage failed for $BAM" >&2

done < "$LIST"
