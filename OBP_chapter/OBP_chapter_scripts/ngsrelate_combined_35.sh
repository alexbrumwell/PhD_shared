#!/bin/bash
#SBATCH --job-name=ngsrelate_combined_35
#SBATCH --error=jobs/%x_%j.err
#SBATCH --output=jobs/%x_%j.log
#SBATCH --time=24:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# ========= MODULES =========
module purge

module load angsd/0.940 ngsrelate/2.0  # or whatever the module is called

# VARIABLES
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
BAM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_35.list"
AUTOSOMES="/projects/erode/people/qhk515/01_OBP/data/OBP_scaffolds_ge100kb.list"
OUTDIR="/projects/erode/scratch/OBP/ngsrelate/"
FREQS="${OUTDIR}/OBP_combined_35.freqs"
RESULTS_PREFIX="${OUTDIR}/OBP_combined_35"
mkdir -p $OUTDIR

# Minimum individuals
N=$(wc -l < "${BAM_LIST}")
PERCENT_IND=90
MIN_IND=$(( (N * PERCENT_IND + 99) / 100 ))

# ========= STEP 1: ANGSD – genotype likelihoods + MAFs =========
echo "$(timestamp) Starting ANGSD for $RESULTS_PREFIX"

### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
angsd -nThreads $SLURM_CPUS_PER_TASK -ref ${REF} -rf $AUTOSOMES -bam $BAM_LIST -out ${RESULTS_PREFIX} -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -C 50  -only_proper_pairs 1 -skipTriallelic 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -rmTrans 1

echo "$(timestamp) Finished running ANGSD"

N_SITES=$(zcat ${RESULTS_PREFIX}.mafs.gz | tail -n +2 | wc -l)
echo "Number of variable sites (no header): ${N_SITES}"


# ========= STEP 2: Make frequency file for NgsRelate =========
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat ${RESULTS_PREFIX}.mafs.gz | cut -f 6 | sed 1d > $FREQS


# ========= STEP 3: Make ID file (sample labels) =========
# To ensure that NgsRelate results have sample labels
IDS_FILE="${RESULTS_PREFIX}.ids"

echo "[$(timestamp)] Creating IDs file from BAM list…"

# If it’s just paths, you can use basename instead.
for line in $(cat "${BAM_LIST}"); do
  basename "${line}" ".rmdup.bam"
done > "${IDS_FILE}"


### run NgsRelate
ngsRelate -g ${RESULTS_PREFIX}.glf.gz -p $SLURM_CPUS_PER_TASK -n $N -f $FREQS -O ${RESULTS_PREFIX}.ngsrelate -z $IDS_FILE

echo "$(timestamp) Finished running NGSrelate"
echo "Done"





exit 0

#Later in R

library(tidyverse)
library(ggplot2)
library("Matrix")
library("corrplot")

setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Relatedness")


list.files() #see what files are available to read

# Read in files
base_name <- "OBP_combined_35"
#base_name <- "OBP_mod_19"
#base_name <- "OBP_hist_16"


rel_ce <- read.table(paste0(base_name, ".ngsrelate"), header = TRUE)
Names <- scan(paste0(base_name,".ids"), what = "character")
N_sites <- rel_ce$nSites[1]

# NGSrelate indices are 0-based; convert to 1-based for R
rel_ce$a <- rel_ce$a + 1L
rel_ce$b <- rel_ce$b + 1L

# Build theta matrix
n <- length(Names)

mat <- sparseMatrix(
  i        = rel_ce$a,
  j        = rel_ce$b,
  x        = rel_ce$theta,
  dims     = c(n, n),
  symmetric = TRUE
)

diag(mat) <- 0.5   # self-relatedness
mat <- as.matrix(mat)

colnames(mat) <- Names
rownames(mat) <- Names


# Plot relatedness heatmap
pdf(paste0(base_name,"_relatedness.pdf"), width = 8, height = 8)
corrplot(
  mat,
  tl.col  = "black",
  type    = "lower",
  tl.cex  = 0.95,
  tl.srt  = 45,
  pch.cex = 1.0,
  method  = "square",
  is.corr = FALSE
)
mtext(
  paste0("Relatedness among ", base_name, " OBP individuals (N_sites = ", N_sites, ")"),
  side = 1,   # bottom
  line = 4,
  cex = 1.1
)
dev.off()


# Extract related samples
rel_ce$name_a <- Names[rel_ce$a]
rel_ce$name_b <- Names[rel_ce$b]

# First-degree or closer (approx. theta >= 0.1875)
first_deg <- subset(rel_ce, theta >= 0.1875)
first_deg[, c("name_a", "name_b", "theta")]

# Second-degree (approx. 0.09375–0.1875)
second_deg <- subset(rel_ce,
                     theta >= 0.09375 & theta < 0.1875)
second_deg[, c("name_a", "name_b", "theta")]




