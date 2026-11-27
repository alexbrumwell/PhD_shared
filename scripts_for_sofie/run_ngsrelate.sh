#!/bin/bash
#SBATCH --job-name=ngsrelate
#SBATCH --error=jobs/%x_%j.err
#SBATCH --output=jobs/%x_%j.log
#SBATCH --time=24:00:00
#SBATCH --partition=cpuqueue                 
#SBATCH --qos=normal                         
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@kent.ac.uk    # EDIT

# ==============================================================================
# PAIRWISE RELATEDNESS WITH ANGSD + NgsRelate
#
# What this script does:
#   STEP 1. ANGSD: genotype likelihoods and allele frequencies at SNPs.
#   STEP 2. Extracts the allele frequency column in the format NgsRelate needs.
#   STEP 3. Makes a sample ID file from the BAM list, so results have names.
#   STEP 4. Runs NgsRelate.
#
# How to run:
#   1. Edit the USER SETTINGS block below (or pass values on the command line).
#   2. Make the log folder first (Slurm will not create it):   mkdir -p jobs
#   3. Submit:
#        sbatch run_ngsrelate.sh
#      or, to run several sample subsets with the same settings:
#        sbatch run_ngsrelate.sh all_samples.list  combined
#        sbatch run_ngsrelate.sh modern.list       modern
#        sbatch run_ngsrelate.sh historical.list   historical
#      Argument 1 = BAM list, argument 2 = name used for the output files.
#
# Inputs:
#   REF        reference genome FASTA (indexed with samtools faidx)
#   BAM_LIST   text file with one BAM path per line
#   REGIONS    text file with one scaffold/chromosome name per line, i.e. the
#              scaffolds to analyse. This should normally be autosomes only:
#              leave out sex-linked scaffolds, the mitogenome and short
#              scaffolds (e.g. < 100 kb).
#
# Main output (in OUTDIR):
#   <name>.ngsrelate   NgsRelate results; one row per pair of samples
#   <name>.ids         sample names, in the same order as the BAM list
#   <name>.freqs       allele frequencies given to NgsRelate
#   Use plot_ngsrelate.R afterwards to plot the results and list related pairs.
# ==============================================================================

set -euo pipefail
timestamp() { date '+%F %T'; }


# ============================================================
# USER SETTINGS (EDIT)
# ============================================================

# Module names differ between clusters. Check with:  module avail angsd
ANGSD_MODULE="angsd/0.940"
NGSRELATE_MODULE="ngsrelate/2.0"

REF="/path/to/reference.fa"                          # EDIT
REGIONS="/path/to/autosomal_scaffolds.list"          # EDIT
OUTDIR="/path/to/output/ngsrelate"                   # EDIT

# Default BAM list and output name. Command-line arguments override these.
BAM_LIST="${1:-/path/to/bam.list}"                   # EDIT
NAME="${2:-combined}"                                # EDIT

# A site is kept only if at least this % of individuals have data there.
PERCENT_IND=90

# SNP and read filters (the values used in the OBP analysis).
MIN_MAPQ=30         # minimum mapping quality
MIN_BASEQ=20        # minimum base quality
MIN_MAF=0.05        # minimum minor allele frequency
SNP_PVAL=1e-6       # SNP p-value cut-off
RM_TRANS=1          # 1 = remove transitions (recommended with historical/
                    #     damaged DNA); 0 = keep them (modern-only data)
PROPER_PAIRS=1      # 1 = only properly paired reads. If your historical reads
                    #     were collapsed (merged) into single reads, check that
                    #     they are not being discarded; set to 0 if so.


# ============================================================
# CHECKS (stop early with a clear message if something is wrong)
# ============================================================
module purge
module load "${ANGSD_MODULE}" "${NGSRELATE_MODULE}"

# Use the Slurm CPU count, or 1 if the script is run outside Slurm.
THREADS="${SLURM_CPUS_PER_TASK:-1}"

for f in "${REF}" "${REGIONS}" "${BAM_LIST}"; do
    [[ -s "${f}" ]] || { echo "ERROR: file missing or empty: ${f}"; exit 1; }
done
[[ -s "${REF}.fai" ]] || { echo "ERROR: reference not indexed. Run: samtools faidx ${REF}"; exit 1; }

mkdir -p "${OUTDIR}"
PREFIX="${OUTDIR%/}/${NAME}"
FREQS="${PREFIX}.freqs"
IDS_FILE="${PREFIX}.ids"

# Clean copy of the BAM list: no blank lines, no Windows line endings.
# (Blank lines would make the sample count wrong.)
CLEAN_LIST="${PREFIX}.bamlist"
tr -d '\r' < "${BAM_LIST}" | awk 'NF' > "${CLEAN_LIST}"

# Check that every BAM exists.
MISSING=0
while IFS= read -r bam; do
    if [[ ! -s "${bam}" ]]; then echo "ERROR: BAM not found: ${bam}"; MISSING=1; fi
done < "${CLEAN_LIST}"
[[ ${MISSING} -eq 0 ]] || exit 1

N=$(wc -l < "${CLEAN_LIST}")
MIN_IND=$(( (N * PERCENT_IND + 99) / 100 ))   # rounded up

echo "$(timestamp) Run name:     ${NAME}"
echo "$(timestamp) Samples:      ${N}"
echo "$(timestamp) Min samples per site (${PERCENT_IND}%): ${MIN_IND}"
echo "$(timestamp) Threads:      ${THREADS}"


# ============================================================
# STEP 1: ANGSD - genotype likelihoods + allele frequencies
# ============================================================
echo "$(timestamp) Starting ANGSD"

angsd \
    -nThreads "${THREADS}" \
    -ref "${REF}" \
    -rf "${REGIONS}" \
    -bam "${CLEAN_LIST}" \
    -out "${PREFIX}" \
    -gl 2 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -snp_pval "${SNP_PVAL}" \
    -minMaf "${MIN_MAF}" \
    -doGlf 3 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs "${PROPER_PAIRS}" \
    -C 50 \
    -skipTriallelic 1 \
    -minMapQ "${MIN_MAPQ}" \
    -minQ "${MIN_BASEQ}" \
    -minInd "${MIN_IND}" \
    -rmTrans "${RM_TRANS}"

echo "$(timestamp) Finished ANGSD"

N_SITES=$(zcat "${PREFIX}.mafs.gz" | tail -n +2 | wc -l)
echo "$(timestamp) Variable sites kept: ${N_SITES}"
if [[ "${N_SITES}" -eq 0 ]]; then
    echo "ERROR: no sites passed the filters. Try lowering PERCENT_IND,"
    echo "  or check that the scaffold names in REGIONS match the reference."
    exit 1
fi


# ============================================================
# STEP 2: Frequency file for NgsRelate
# ============================================================
# Find the allele frequency column by its NAME (knownEM) rather than its
# number, because the column number changes depending on ANGSD options.
zcat "${PREFIX}.mafs.gz" | awk '
NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == "knownEM") col = i
    if (!col) { print "ERROR: no knownEM column in mafs file" > "/dev/stderr"; exit 1 }
    next
}
{ print $col }
' > "${FREQS}"


# ============================================================
# STEP 3: Sample ID file
# ============================================================
# Sample name = BAM file name without folders and without everything from
# the first "." onwards, e.g.  /data/OBP_01.rmdup.sorted.bam  ->  OBP_01
# Edit the sed command if your file names need something different.
while IFS= read -r bam; do
    basename "${bam}" | sed 's/\..*$//'
done < "${CLEAN_LIST}" > "${IDS_FILE}"

# Warn if two BAMs give the same sample name.
DUPS=$(sort "${IDS_FILE}" | uniq -d)
if [[ -n "${DUPS}" ]]; then
    echo "WARNING: duplicate sample names in ${IDS_FILE}:"
    echo "${DUPS}"
fi


# ============================================================
# STEP 4: NgsRelate
# ============================================================
echo "$(timestamp) Starting NgsRelate"

ngsRelate \
    -g "${PREFIX}.glf.gz" \
    -p "${THREADS}" \
    -n "${N}" \
    -f "${FREQS}" \
    -z "${IDS_FILE}" \
    -O "${PREFIX}.ngsrelate"

echo "$(timestamp) Finished NgsRelate"
echo "Results: ${PREFIX}.ngsrelate"
echo "Sample IDs: ${IDS_FILE}"
echo "Next: copy both files to your computer and run plot_ngsrelate.R"
exit 0



