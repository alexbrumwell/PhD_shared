#!/bin/bash
#SBATCH --job-name=sex_scaffolds
#SBATCH --partition=cpuqueue                 
#SBATCH --qos=normal    
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@kent.ac.uk    # EDIT

# NOTE: Slurm will not create the jobs/ folder for you. Run "mkdir -p jobs"
# in the folder you submit from BEFORE running "sbatch".

# SCRIPT TO DETERMINE SEX AND IDENTIFY SEX-LINKED SCAFFOLDS
# For any bird (ZW) whole-genome dataset, including temporal datasets that
# mix historical and modern samples. Afterwards, plot the results with
# plot_sex_scaffolds.R.

# ==============================================================================
# This script:
#
#   STAGE 1. Runs samtools coverage on every BAM file (mean depth per scaffold).
#            Existing coverage files are re-used, so re-running is fast.
#
#   STAGE 2. Finds the Z-linked scaffold used for sexing.
#            - If SEX_SCAFFOLD is set below, that scaffold is used.
#            - If SEX_SCAFFOLD is left empty (""), the script searches for it:
#              it looks for large scaffolds where samples split into two
#              groups, one near 0.5x autosomal depth (ZW females) and one near
#              1.0x (ZZ males). An autosome has every sample near 1.0, so it
#              does not show this split. The largest scaffold that passes is
#              used. ALWAYS check this choice (see z_scaffold_search.tsv).
#
#   STAGE 3. For each sample:
#            - calculates a length-weighted mean autosomal depth,
#            - normalises the Z scaffold depth against it,
#            - assigns a provisional sex:
#                  male_ZZ:       normalised Z depth near 1.0
#                  female_ZW:     normalised Z depth near 0.5
#                  ambiguous:     normalised Z depth below the female range
#                  low_coverage:  autosomal depth below MIN_AUTO_DEPTH
#            - normalises the depth of every scaffold.
#
#   STAGE 4. Compares mean scaffold depth between provisional males and
#            females to find additional candidate Z- and W-linked scaffolds.
#            Ambiguous and low_coverage samples are left out of this step.
#
# Input:
#   - All BAM files located directly within BAM_DIR.
#
# Main outputs (in OUT_DIR):
#   z_scaffold_search.tsv            (only if SEX_SCAFFOLD was left empty)
#   sex_results.tsv
#   all_normalised_depths.tsv
#   possible_sex_linked_scaffolds.tsv
#
# Suggested workflow:
#   Run 1: leave SEX_SCAFFOLD="" and EXCLUDE_LIST="".
#   Check the outputs. Put all Z/W candidates you accept into a text file
#   (one scaffold name per line) and set EXCLUDE_LIST to it.
#   Run 2: re-run. Those scaffolds are then removed from the autosomal
#   baseline, which makes the female ~0.5 / male ~1.0 values cleaner.
# ==============================================================================

set -euo pipefail
shopt -s nullglob

export LC_NUMERIC=C # Standardise numeric formatting between Bash and awk.

module load samtools/1.21


# ============================================================
# USER SETTINGS (EDIT)
# ============================================================

# Directory containing the input BAM files (historical AND modern together).
BAM_DIR="/path/to/bams"                                              # EDIT

# Directory in which all results and intermediate files will be written.
OUT_DIR="/path/to/sex_scaffold_analysis"                             # EDIT

# Known Z-linked scaffold used for the initial sex assignment.
# Leave empty ("") to let the script search for it (STAGE 2).
# If the reference is chromosome-level, the Z may already be named
# (check with: samtools view -H your.bam | grep "^@SQ" | head -50).
SEX_SCAFFOLD=""

# Mitochondrial scaffold, excluded from the autosomal baseline.
# Optional: leave empty if not known. The mitogenome is ~17 kb, so it is
# already removed by MIN_LENGTH unless it is part of a longer scaffold.
MITO_SCAFFOLD=""

# Optional: file with one scaffold name per line to ALSO exclude from the
# autosomal baseline (e.g. Z/W candidates accepted after the first run).
EXCLUDE_LIST=""

# Minimum scaffold length (bp) used for the autosomal baseline and the
# candidate sex-linked scaffold screen.
MIN_LENGTH=100000

# Minimum scaffold length (bp) for the automatic Z scaffold search.
# Larger scaffolds give more reliable depth estimates.
DETECT_MIN_LENGTH=1000000

# Minimum mapping quality passed to samtools coverage.
# 0 keeps every read. A value of 20-30 removes multi-mapping reads, which
# matters for Z/W gametolog regions and short historical reads.
# Use the same value for every sample.
MIN_MAPQ=20

# Samples with autosomal depth below this value are flagged "low_coverage"
# and left out of the male/female comparison (STAGE 4).
MIN_AUTO_DEPTH=1

# Sex thresholds for normalised Z depth (Z depth / autosomal depth).
MALE_MIN=0.75      # >= this value: male_ZZ
FEMALE_MIN=0.30    # >= this value and < MALE_MIN: female_ZW; below: ambiguous


# ============================================================
# PREPARE OUTPUT
# ============================================================
mkdir -p "${OUT_DIR}/coverage" "${OUT_DIR}/normalised"

SEX_RESULTS="${OUT_DIR}/sex_results.tsv"
ALL_DEPTHS="${OUT_DIR}/all_normalised_depths.tsv"
Z_SEARCH="${OUT_DIR}/z_scaffold_search.tsv"
SEX_LINKED_RESULTS="${OUT_DIR}/possible_sex_linked_scaffolds.tsv"

# Collect the BAM files and stop early if there are none.
BAMS=( "${BAM_DIR}"/*.bam )
if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: no .bam files found in ${BAM_DIR}"
    exit 1
fi
echo "Found ${#BAMS[@]} BAM files"

# Build one file listing every scaffold to leave out of the autosomal
# baseline (mito + user exclude list). The Z scaffold is added in STAGE 3.
BASE_EXCLUDE="${OUT_DIR}/baseline_exclude.txt"
: > "${BASE_EXCLUDE}"
[[ -n "${MITO_SCAFFOLD}" ]] && echo "${MITO_SCAFFOLD}" >> "${BASE_EXCLUDE}"
if [[ -n "${EXCLUDE_LIST}" ]]; then
    if [[ ! -f "${EXCLUDE_LIST}" ]]; then
        echo "ERROR: EXCLUDE_LIST file not found: ${EXCLUDE_LIST}"
        exit 1
    fi
    # Strip Windows line endings and blank lines in case the file was
    # written on a laptop.
    tr -d '\r' < "${EXCLUDE_LIST}" | awk 'NF' >> "${BASE_EXCLUDE}"
fi


# ------------------------------------------------------------
# Helper: length-weighted autosomal depth for one coverage file.
#   depth = sum(mean_depth x length) / sum(length)
# over scaffolds >= MIN_LENGTH that are NOT in the exclude file.
# ------------------------------------------------------------
autosomal_depth () {
    local coverage_file=$1
    local exclude_file=$2
    awk -v min_length="${MIN_LENGTH}" '
    BEGIN { FS = "\t" }
    # First file: the scaffolds to exclude. (Checked by file name rather
    # than NR == FNR, because that trick breaks when the file is empty.)
    FILENAME == ARGV[1] { if (NF) exclude[$1] = 1; next }
    # Skip the samtools coverage header.
    /^#/ { next }
    {
        scaffold        = $1
        scaffold_length = $3 - $2 + 1
        mean_depth      = $7
        if (!(scaffold in exclude) && scaffold_length >= min_length) {
            depth_sum  += mean_depth * scaffold_length
            length_sum += scaffold_length
        }
    }
    END {
        if (length_sum > 0 && depth_sum > 0) printf "%.6f", depth_sum / length_sum
        else print "NA"
    }' "${exclude_file}" "${coverage_file}"
}


# ============================================================
# STAGE 1: COVERAGE FOR EACH BAM
# ============================================================
for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "${BAM}" .bam)
    COVERAGE_FILE="${OUT_DIR}/coverage/${SAMPLE}.coverage.tsv"

    # Re-use an existing coverage file. Delete the coverage/ folder if you
    # change MIN_MAPQ, so everything is recalculated with the new value.
    if [[ -s "${COVERAGE_FILE}" ]]; then
        echo "Coverage exists for ${SAMPLE}, skipping samtools"
        continue
    fi

    echo "Running samtools coverage for ${SAMPLE}"

    if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" && ! -f "${BAM}.csi" ]]; then
        echo "  Indexing ${SAMPLE}"
        samtools index "${BAM}"
    fi

    # Write to a temporary file first so a failed run does not leave a
    # half-written file that would be "re-used" next time.
    samtools coverage -q "${MIN_MAPQ}" "${BAM}" > "${COVERAGE_FILE}.tmp"
    mv "${COVERAGE_FILE}.tmp" "${COVERAGE_FILE}"
done


# ============================================================
# STAGE 2: CHOOSE THE Z SCAFFOLD
# ============================================================
if [[ -z "${SEX_SCAFFOLD}" ]]; then
    echo
    echo "SEX_SCAFFOLD not set: searching for the Z scaffold"

    # Step 2a: a provisional normalised depth for every large scaffold in
    # every sample (baseline = all long scaffolds, minus mito/exclude list).
    PRELIM="${OUT_DIR}/prelim_normalised.tsv"
    : > "${PRELIM}"
    for BAM in "${BAMS[@]}"; do
        SAMPLE=$(basename "${BAM}" .bam)
        COVERAGE_FILE="${OUT_DIR}/coverage/${SAMPLE}.coverage.tsv"
        AUTO=$(autosomal_depth "${COVERAGE_FILE}" "${BASE_EXCLUDE}")
        if [[ "${AUTO}" == "NA" ]]; then
            echo "  WARNING: no usable depth for ${SAMPLE}, left out of search"
            continue
        fi
        # Very low coverage samples give noisy ratios, so leave them out.
        if awk -v d="${AUTO}" -v m="${MIN_AUTO_DEPTH}" 'BEGIN{exit !(d < m)}'; then
            echo "  ${SAMPLE}: autosomal depth ${AUTO} < ${MIN_AUTO_DEPTH}, left out of search"
            continue
        fi
        awk -v sample="${SAMPLE}" -v auto="${AUTO}" -v min_len="${DETECT_MIN_LENGTH}" '
        BEGIN { FS = OFS = "\t" }
        /^#/ { next }
        ($3 - $2 + 1) >= min_len { print sample, $1, $3 - $2 + 1, $7 / auto }
        ' "${COVERAGE_FILE}" >> "${PRELIM}"
    done

    # Step 2b: for each scaffold, count samples in the "female" band
    # (0.35-0.70) and the "male" band (0.80-1.25).
    #   Z scaffold: most samples fall in one of the two bands, and BOTH
    #               bands contain at least 2 samples.
    #   Autosome:   almost all samples in the male band only -> rejected.
    #   W scaffold: females ~0.5, males ~0 -> few in male band -> rejected.
    printf "Scaffold\tLength\tN_samples\tN_female_band\tN_male_band\tFraction_in_bands\tPasses\n" \
        > "${Z_SEARCH}"
    awk '
    BEGIN { FS = OFS = "\t" }
    {
        scaffold = $2
        len[scaffold] = $3
        n[scaffold]++
        if ($4 >= 0.35 && $4 <= 0.70) low[scaffold]++
        else if ($4 >= 0.80 && $4 <= 1.25) high[scaffold]++
    }
    END {
        for (s in n) {
            l = low[s] + 0; h = high[s] + 0
            frac = (l + h) / n[s]
            pass = (l >= 2 && h >= 2 && frac >= 0.80) ? "yes" : "no"
            printf "%s\t%s\t%d\t%d\t%d\t%.3f\t%s\n", s, len[s], n[s], l, h, frac, pass
        }
    }' "${PRELIM}" | sort -t $'\t' -k7,7r -k2,2nr >> "${Z_SEARCH}"

    # Step 2c: pick the longest scaffold that passes.
    SEX_SCAFFOLD=$(awk -F'\t' 'NR > 1 && $7 == "yes" { print $1; exit }' "${Z_SEARCH}")

    if [[ -z "${SEX_SCAFFOLD}" ]]; then
        echo "ERROR: no scaffold showed a clear male/female depth split."
        echo "  Possible reasons: all samples are one sex, coverage is too low,"
        echo "  or the Z is split into scaffolds < ${DETECT_MIN_LENGTH} bp."
        echo "  Check ${Z_SEARCH}, or set SEX_SCAFFOLD by hand."
        exit 1
    fi

    N_PASS=$(awk -F'\t' 'NR > 1 && $7 == "yes"' "${Z_SEARCH}" | wc -l)
    echo "  ${N_PASS} scaffold(s) passed; using the longest: ${SEX_SCAFFOLD}"
    echo "  Check ${Z_SEARCH} before trusting this choice."
else
    echo "Using user-specified sex scaffold: ${SEX_SCAFFOLD}"
fi


# ============================================================
# STAGE 3: SEX ASSIGNMENT AND NORMALISED DEPTH PER SAMPLE
# ============================================================

# Exclude list for the final baseline = mito + user list + Z scaffold.
FINAL_EXCLUDE="${OUT_DIR}/final_baseline_exclude.txt"
cat "${BASE_EXCLUDE}" > "${FINAL_EXCLUDE}"
echo "${SEX_SCAFFOLD}" >> "${FINAL_EXCLUDE}"

printf "Sample\tAutosomal_depth\tSex_scaffold\tSex_scaffold_depth\tNormalised_depth\tPredicted_sex\n" \
    > "${SEX_RESULTS}"
printf "Sample\tScaffold\tLength\tMean_depth\tAutosomal_depth\tNormalised_depth\n" \
    > "${ALL_DEPTHS}"

for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "${BAM}" .bam)
    COVERAGE_FILE="${OUT_DIR}/coverage/${SAMPLE}.coverage.tsv"
    NORMALISED_FILE="${OUT_DIR}/normalised/${SAMPLE}.normalised.tsv"

    echo
    echo "Processing ${SAMPLE}"

    AUTOSOMAL_DEPTH=$(autosomal_depth "${COVERAGE_FILE}" "${FINAL_EXCLUDE}")

    # One bad sample is skipped rather than stopping the whole run.
    if [[ "${AUTOSOMAL_DEPTH}" == "NA" ]]; then
        echo "  WARNING: could not calculate autosomal depth, skipping ${SAMPLE}"
        continue
    fi

    SEX_DEPTH=$(awk -v s="${SEX_SCAFFOLD}" -F'\t' '!/^#/ && $1 == s { print $7 }' "${COVERAGE_FILE}")
    if [[ -z "${SEX_DEPTH}" ]]; then
        echo "  WARNING: ${SEX_SCAFFOLD} not found, skipping ${SAMPLE}"
        continue
    fi

    NORMALISED_SEX_DEPTH=$(awk -v z="${SEX_DEPTH}" -v a="${AUTOSOMAL_DEPTH}" \
        'BEGIN { printf "%.4f", z / a }')

    # Provisional sex. low_coverage is checked first: at very low depth the
    # ratio is too noisy to trust.
    PREDICTED_SEX=$(awk \
        -v ratio="${NORMALISED_SEX_DEPTH}" -v auto="${AUTOSOMAL_DEPTH}" \
        -v min_auto="${MIN_AUTO_DEPTH}" \
        -v male_min="${MALE_MIN}" -v female_min="${FEMALE_MIN}" '
    BEGIN {
        if (auto < min_auto)        print "low_coverage"
        else if (ratio >= male_min)   print "male_ZZ"
        else if (ratio >= female_min) print "female_ZW"
        else                          print "ambiguous"
    }')

    printf "%s\t%.4f\t%s\t%.4f\t%.4f\t%s\n" \
        "${SAMPLE}" "${AUTOSOMAL_DEPTH}" "${SEX_SCAFFOLD}" \
        "${SEX_DEPTH}" "${NORMALISED_SEX_DEPTH}" "${PREDICTED_SEX}" \
        >> "${SEX_RESULTS}"

    # Normalised depth of every scaffold for this sample.
    printf "Sample\tScaffold\tLength\tMean_depth\tAutosomal_depth\tNormalised_depth\n" \
        > "${NORMALISED_FILE}"
    awk -v sample="${SAMPLE}" -v auto="${AUTOSOMAL_DEPTH}" '
    BEGIN { FS = OFS = "\t" }
    /^#/ { next }
    { print sample, $1, $3 - $2 + 1, $7, auto, $7 / auto }
    ' "${COVERAGE_FILE}" >> "${NORMALISED_FILE}"

    tail -n +2 "${NORMALISED_FILE}" >> "${ALL_DEPTHS}"

    echo "  Autosomal depth:  ${AUTOSOMAL_DEPTH}"
    echo "  ${SEX_SCAFFOLD} depth: ${SEX_DEPTH}"
    echo "  Normalised depth: ${NORMALISED_SEX_DEPTH}"
    echo "  Predicted sex:    ${PREDICTED_SEX}"
done

echo
echo "Sex results: ${SEX_RESULTS}"
echo "All normalised scaffold depths: ${ALL_DEPTHS}"

# Summary of how many samples ended up in each class.
echo "Sex class counts:"
awk -F'\t' 'NR > 1 { c[$6]++ } END { for (k in c) printf "  %s\t%d\n", k, c[k] }' "${SEX_RESULTS}"


# ============================================================
# STAGE 4: SEARCH FOR OTHER POSSIBLE SEX-LINKED SCAFFOLDS
# ============================================================

N_MALE=$(awk -F'\t' 'NR > 1 && $6 == "male_ZZ"' "${SEX_RESULTS}" | wc -l)
N_FEMALE=$(awk -F'\t' 'NR > 1 && $6 == "female_ZW"' "${SEX_RESULTS}" | wc -l)
if [[ ${N_MALE} -eq 0 || ${N_FEMALE} -eq 0 ]]; then
    echo "WARNING: need at least one male and one female for the scaffold screen"
    echo "  (males: ${N_MALE}, females: ${N_FEMALE}). Skipping STAGE 4."
    exit 0
fi

awk -v min_length="${MIN_LENGTH}" '
BEGIN { FS = OFS = "\t" }

# First file: predicted sex per sample.
NR == FNR {
    if (FNR > 1) sample_sex[$1] = $6
    next
}

# Second file: skip its header.
FNR == 1 { next }

{
    sample = $1; scaffold = $2
    scaffold_length[scaffold] = $3
    normalised_depth = $6

    # ambiguous and low_coverage samples are not used.
    if (sample_sex[sample] == "male_ZZ") {
        male_sum[scaffold] += normalised_depth; male_count[scaffold]++
    } else if (sample_sex[sample] == "female_ZW") {
        female_sum[scaffold] += normalised_depth; female_count[scaffold]++
    }
}

END {
    print "Scaffold", "Length", "N_male", "N_female", "Male_mean", \
          "Female_mean", "Female_male_ratio", "Classification"

    for (scaffold in scaffold_length) {
        if (scaffold_length[scaffold] < min_length) continue
        if (male_count[scaffold] == 0 || female_count[scaffold] == 0) continue

        male_mean   = male_sum[scaffold] / male_count[scaffold]
        female_mean = female_sum[scaffold] / female_count[scaffold]

        ratio = "NA"
        classification = "autosomal_or_unknown"
        if (male_mean > 0) ratio = female_mean / male_mean

        # Z-linked: males ~1, females ~0.5, ratio ~0.5.
        if (ratio != "NA" &&
            male_mean >= 0.75 && male_mean <= 1.25 &&
            female_mean >= 0.25 && female_mean <= 0.75 &&
            ratio >= 0.30 && ratio <= 0.70) {
            classification = "possible_Z_linked"
        }
        # W-linked: little or no depth in males, some in females.
        else if (male_mean <= 0.10 && female_mean >= 0.20) {
            classification = "possible_W_linked"
        }

        if (ratio != "NA") ratio = sprintf("%.4f", ratio)
        printf "%s\t%s\t%d\t%d\t%.4f\t%.4f\t%s\t%s\n", \
               scaffold, scaffold_length[scaffold], \
               male_count[scaffold], female_count[scaffold], \
               male_mean, female_mean, ratio, classification
    }
}
' "${SEX_RESULTS}" "${ALL_DEPTHS}" |
{
    # Keep the header, then sort by classification and length (largest first).
    IFS= read -r header
    printf '%s\n' "${header}"
    sort -t $'\t' -k8,8 -k2,2nr
} > "${SEX_LINKED_RESULTS}"

echo
echo "Possible sex-linked scaffolds: ${SEX_LINKED_RESULTS}"
echo "Candidate counts:"
awk -F'\t' 'NR > 1 && $8 != "autosomal_or_unknown" { c[$8]++ }
            END { for (k in c) printf "  %s\t%d\n", k, c[k] }' "${SEX_LINKED_RESULTS}"

exit 0
