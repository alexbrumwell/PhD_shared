#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal 
#SBATCH --cpus-per-task=4
#SBATCH --job-name=OBPs_coverage
#SBATCH --output=/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/scripts/OBPs_samtoolscov_%j.out
#SBATCH --error=/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/scripts/OBPs_samtoolscov_%j.err
#SBATCH --time=24:00:00

# No need to load module as qualimap is downloaded in my directory
module load samtools


# Set variables
DIR="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/"


for file in ${DIR}*.bam; do
    samtools coverage ${file} > ${DIR}${file%.bam}_coverage.tsv
    echo "Coverage for ${file}:"  >> ${DIR}OBP_bam_coverage_histogram.txt
    samtools coverage ${file} --histogram >> ${DIR}OBP_bam_coverage_histogram.txt
    echo "Depth for ${file}:"  >> ${DIR}OBP_bam_coverage_depth.txt
    samtools coverage ${file} --plot-depth >> ${DIR}OBP_bam_coverage_depth.txt
done

