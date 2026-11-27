#!/bin/bash
#SBATCH --job-name=sex_scaffolds
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

# SCRIPT TO DETERMINE SEX AND IDENTIFY SEX-LINKED SCAFFOLDS

# ==============================================================================
# This script:
#
#   1. Calculates mean sequencing depth for every scaffold in each BAM file.
#   2. Calculates a length-weighted mean autosomal depth for each sample.
#   3. Normalises the depth of a known or suspected sex-linked scaffold (in this study scaffold_5)
#      against autosomal depth.
#   4. Assigns provisional avian sex classes:
#          male_ZZ:    normalised Z depth expected near 1.0
#          female_ZW:  normalised Z depth expected near 0.5
#          ambiguous:  normalised Z depth below the expected female range
#   5. Normalises the depth of every scaffold within each sample.
#   6. Compares mean scaffold depth between provisional males and females to
#      identify additional candidate Z- and W-linked scaffolds.

# Input:
#   - All BAM files located directly within BAM_DIR.

# Main outputs:
#
#   sex_results.tsv
#       Per-sample autosomal depth, sex-scaffold depth, normalised depth,
#       and provisional sex assignment.
#
#   all_normalised_depths.tsv
#       Mean and normalised depth for every scaffold in every sample.
#
#   possible_sex_linked_scaffolds.tsv
#       Mean normalised scaffold depth in provisional males and females,
#       the female-to-male depth ratio, and scaffold classification.
# ==============================================================================

set -euo pipefail

export LC_NUMERIC=C # Standardise numeric formatting between Bash and awk.

module load samtools/1.21


# ============================================================
# USER SETTINGS (EDIT)
# ============================================================

# Directory containing the input BAM files.
BAM_DIR="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP"

# Directory in which all results and intermediate files will be written.
OUT_DIR="/projects/erode/scratch/OBP/sex_scaffold_analysis"

# Known or suspected Z-linked scaffold used for the initial sex assignment.
SEX_SCAFFOLD="scaffold_5"

# Mitochondrial scaffold, excluded from the autosomal depth calculation
MITO_SCAFFOLD="scaffold_1904"

# Minimum scaffold length, in base pairs
MIN_LENGTH=100000

# Minimum mapping quality passed to samtools coverage.
MIN_MAPQ=0

# Female if scaffold_5/autosomal depth is below this value
THRESHOLD=0.75


# ============================================================
# PREPARE OUTPUT
# ============================================================
mkdir -p "${OUT_DIR}/coverage"
mkdir -p "${OUT_DIR}/normalised"

# Main table containing the provisional sex assignment for every sample.
SEX_RESULTS="${OUT_DIR}/sex_results.tsv"

# Combined table containing the normalised depth of every scaffold in every sample.
ALL_DEPTHS="${OUT_DIR}/all_normalised_depths.tsv"

# Initialise the sex-results table and write its column names.
#
# Columns:
#   Sample            Sample identifier derived from the BAM filename
#   Autosomal_depth   Length-weighted mean depth across selected autosomes
#   Scaffold_5_depth  Mean depth of the putative Z-linked scaffold
#   Normalised_depth  Scaffold_5 depth divided by autosomal depth
#   Predicted_sex     Provisional sex assignment
printf "Sample\tAutosomal_depth\tScaffold_5_depth\tNormalised_depth\tPredicted_sex\n" \
    > "${SEX_RESULTS}"

# Initialise the combined scaffold-depth table and write its column names.
#
# Columns:
#   Sample             Sample identifier
#   Scaffold           Scaffold name
#   Length             Scaffold length in base pairs
#   Mean_depth         Mean sequencing depth reported by samtools coverage
#   Autosomal_depth    Sample-specific autosomal depth
#   Normalised_depth   Scaffold depth divided by autosomal depth
printf "Sample\tScaffold\tLength\tMean_depth\tAutosomal_depth\tNormalised_depth\n" \
    > "${ALL_DEPTHS}"


# ============================================================
# PROCESS EACH BAM
# ============================================================

# Process every BAM file found directly within BAM_DIR.
for BAM in "${BAM_DIR}"/*.bam; do

    # Derive the sample identifier by removing the .bam suffix.
    SAMPLE=$(basename "${BAM}" .bam)

    echo "Processing ${SAMPLE}"

    # Per-sample samtools coverage output.
    COVERAGE_FILE="${OUT_DIR}/coverage/${SAMPLE}.coverage.tsv"

    # Per-sample normalised scaffold-depth output.
    NORMALISED_FILE="${OUT_DIR}/normalised/${SAMPLE}.normalised.tsv"


    # --------------------------------------------------------
    # Check that the BAM is indexed
    # --------------------------------------------------------
    if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
        echo "Indexing ${SAMPLE}"
        samtools index "${BAM}"
    fi


    # --------------------------------------------------------
    # Calculate mean depth for every scaffold
    # --------------------------------------------------------
    samtools coverage \
        -q "${MIN_MAPQ}" \
        "${BAM}" \
        > "${COVERAGE_FILE}"


    # --------------------------------------------------------
    # Calculate autosomal depth
    # --------------------------------------------------------

    # Calculate a sample-specific autosomal depth baseline.
    # For each eligible scaffold:
    #   depth contribution = scaffold mean depth × scaffold length
    # The final autosomal depth is:
    #   sum(depth × length) / sum(length)
    # This weights long scaffolds more heavily than short scaffolds and
    # approximates mean depth across the selected autosomal sequence.
    AUTOSOMAL_DEPTH=$(
    awk \
        -v sex_scaffold="${SEX_SCAFFOLD}" \
        -v mito_scaffold="${MITO_SCAFFOLD}" \
        -v min_length="${MIN_LENGTH}" '
    BEGIN {
        # samtools coverage output is tab-delimited.
        FS = "\t"
    }

    # Skip the samtools coverage header.
    /^#/ {
        next
    }

    {
        # Extract scaffold name, covered interval length and mean depth.
        scaffold = $1
        scaffold_length = $3 - $2 + 1
        mean_depth = $7

        # Include only sufficiently long presumed autosomal scaffolds.
        if (scaffold != sex_scaffold &&
            scaffold != mito_scaffold &&
            scaffold_length >= min_length) {

            depth_sum += mean_depth * scaffold_length
            length_sum += scaffold_length
        }
    }

    END {
        # Report the length-weighted mean if at least one eligible
        # scaffold was identified.
        if (length_sum > 0) {
            printf "%.6f", depth_sum / length_sum
        } else {
            print "NA"
        }
    }
    ' "${COVERAGE_FILE}"
)

    # Stop the script if no eligible autosomal scaffolds were available.
    if [[ "${AUTOSOMAL_DEPTH}" == "NA" ]]; then
        echo "ERROR: Could not calculate autosomal depth for ${SAMPLE}"
        exit 1
    fi


    # --------------------------------------------------------
    # Get scaffold_5 depth
    # --------------------------------------------------------

    # Extract the mean depth reported in column 7 for the putative
    # Z-linked scaffold.
    SEX_DEPTH=$(
        awk -v sex_scaffold="${SEX_SCAFFOLD}" '
        BEGIN {
            FS="\t"
        }

        !/^#/ && $1 == sex_scaffold {
            print $7
        }
        ' "${COVERAGE_FILE}"
    )


    # Skip the sample if the specified sex scaffold is absent from the
    # samtools coverage output.
    if [[ -z "${SEX_DEPTH}" ]]; then
        echo "WARNING: ${SEX_SCAFFOLD} was not found for ${SAMPLE}"
        continue
    fi


    # --------------------------------------------------------
    # Normalise scaffold_5 depth
    # --------------------------------------------------------

    # Calculate:
    #
    #   normalised Z depth = scaffold_5 mean depth / autosomal mean depth
    #
    # Under the expected avian ZW system:
    #
    #   male ZZ   approximately 1.0
    #   female ZW approximately 0.5
    NORMALISED_SEX_DEPTH=$(
        awk \
            -v sex_depth="${SEX_DEPTH}" \
            -v autosomal_depth="${AUTOSOMAL_DEPTH}" '
        BEGIN {
            printf "%.4f", sex_depth/autosomal_depth
        }
        '
    )


    # --------------------------------------------------------
    # Predict sex
    # --------------------------------------------------------

    # Assign provisional sex using the normalised Z-depth ratio:
    #
    #   ratio >= 0.75:
    #       male_ZZ
    #
    #   ratio >= 0.30 and < 0.75:
    #       female_ZW
    #
    #   ratio < 0.30:
    #       ambiguous
    #
    # The ambiguous category prevents samples with unexpectedly low
    # scaffold_5 coverage from automatically being classified as female.
    PREDICTED_SEX=$(
    awk \
        -v ratio="${NORMALISED_SEX_DEPTH}" '
    BEGIN {
        if (ratio >= 0.75) {
            print "male_ZZ"
        }
        else if (ratio >= 0.30 && ratio < 0.75) {
            print "female_ZW"
        }
        else {
            print "ambiguous"
        }
    }
    '
)


    # --------------------------------------------------------
    # Write sex result
    # --------------------------------------------------------

    # Append the current sample's results to the combined sex-results table.
    printf "%s\t%.4f\t%.4f\t%.4f\t%s\n" \
        "${SAMPLE}" \
        "${AUTOSOMAL_DEPTH}" \
        "${SEX_DEPTH}" \
        "${NORMALISED_SEX_DEPTH}" \
        "${PREDICTED_SEX}" \
        >> "${SEX_RESULTS}"


    # --------------------------------------------------------
    # Normalise every scaffold
    # --------------------------------------------------------

    # Initialise a separate normalised-depth table for the current sample.
    printf "Sample\tScaffold\tLength\tMean_depth\tAutosomal_depth\tNormalised_depth\n" \
        > "${NORMALISED_FILE}"

    # For every scaffold, calculate:
    #
    #   normalised depth = scaffold mean depth / autosomal mean depth
    #
    # This table is subsequently used to compare scaffold-level depth
    # between samples classified as male and female.
    awk -v sample="$SAMPLE" -v autosomal_depth="$AUTOSOMAL_DEPTH" \
'BEGIN{FS=OFS="\t"} /^#/{next} {scaffold=$1; scaffold_length=$3-$2+1; mean_depth=$7; normalised_depth=mean_depth/autosomal_depth; print sample,scaffold,scaffold_length,mean_depth,autosomal_depth,normalised_depth}' \
"$COVERAGE_FILE" >> "$NORMALISED_FILE"

    # Add the sample's scaffold results to the combined table while
    # excluding the per-sample header line.
    tail -n +2 "${NORMALISED_FILE}" >> "${ALL_DEPTHS}"

    # Print a concise summary to the Slurm standard-output log.
    echo "  Autosomal depth: ${AUTOSOMAL_DEPTH}"
    echo "  ${SEX_SCAFFOLD} depth: ${SEX_DEPTH}"
    echo "  Normalised depth: ${NORMALISED_SEX_DEPTH}"
    echo "  Predicted sex: ${PREDICTED_SEX}"

done


# ============================================================
# FINISHED
# ============================================================

# Report the locations of the two main files generated during the
# per-sample processing stage.
echo
echo "Finished."
echo "Sex results: ${SEX_RESULTS}"
echo "All normalised scaffold depths: ${ALL_DEPTHS}"


# =================================================
# SEARCH FOR OTHER POSSIBLE SEX-LINKED SCAFFOLDS
# ============================================================

# Output table containing group-level depth summaries and provisional
# sex-linkage classifications.
SEX_LINKED_RESULTS="${OUT_DIR}/possible_sex_linked_scaffolds.tsv"

# Initialise the candidate-scaffold output.
#
# This file is subsequently overwritten by the awk pipeline below,
# which also writes its own header.
printf "Scaffold\tLength\tMale_mean\tFemale_mean\tFemale_to_male\tClassification\n" \
    > "${SEX_LINKED_RESULTS}"

# Read two input files:
#
#   1. SEX_RESULTS:
#      provides the provisional sex assigned to each sample.
#
#   2. ALL_DEPTHS:
#      provides normalised depth for every sample and scaffold.
#
# The script then calculates mean normalised depth separately for
# provisional males and females.
awk '
BEGIN {
    FS = OFS = "\t"
}

# Read predicted sex from the first file
#
# NR == FNR is true while awk is reading the first input file.
NR == FNR {
    # Skip the sex-results header.
    if (FNR == 1) {
        next
    }

    # Store the predicted sex using the sample identifier as the key.
    sample_sex[$1] = $5
    next
}

# Skip the header of the second file
FNR == 1 {
    next
}

{
    # Extract sample, scaffold, scaffold length and normalised depth.
    sample = $1
    scaffold = $2
    scaffold_length[scaffold] = $3
    normalised_depth = $6

    # Accumulate normalised depth separately for provisional males
    # and provisional females.
    #
    # Samples classified as ambiguous are not included in either group.
    if (sample_sex[sample] == "male_ZZ") {
        male_sum[scaffold] += normalised_depth
        male_count[scaffold]++
    }
    else if (sample_sex[sample] == "female_ZW") {
        female_sum[scaffold] += normalised_depth
        female_count[scaffold]++
    }
}

END {
    # Write the candidate-scaffold table header.
    print "Scaffold", "Length", "Male_mean", "Female_mean", \
          "Female_male_ratio", "Classification"

    # Summarise each scaffold observed in the combined depth table.
    for (scaffold in scaffold_length) {

                # Exclude scaffolds shorter than 100 kb from the
                # candidate sex-linkage screen.
				if (scaffold_length[scaffold] < 100000) {
    			continue
				}

        # A male-to-female comparison cannot be made unless the scaffold
        # has observations from at least one sample in each group.
        if (male_count[scaffold] == 0 || female_count[scaffold] == 0) {
            continue
        }

        # Calculate group mean normalised depth.
        male_mean = male_sum[scaffold] / male_count[scaffold]
        female_mean = female_sum[scaffold] / female_count[scaffold]

        # Default values before applying sex-linkage criteria.
        ratio = "NA"
        classification = "autosomal_or_unknown"

        # Avoid division by zero when calculating the female-to-male ratio.
        if (male_mean > 0) {
            ratio = female_mean / male_mean
        }

        # Candidate Z-linked scaffolds are expected to have:
        #
        #   male normalised depth:   0.75–1.25
        #   female normalised depth: 0.25–0.75
        #   female/male ratio:       0.30–0.70
        if (male_mean >= 0.75 && male_mean <= 1.25 && \
            female_mean >= 0.25 && female_mean <= 0.75 && \
            ratio >= 0.30 && ratio <= 0.70) {

            classification = "possible_Z_linked"
        }

        # Candidate W-linked scaffolds are expected to have:
        #
        #   male normalised depth:   <= 0.10
        #   female normalised depth: >= 0.20
        #
        # This is an exploratory classification and does not itself
        # confirm that a scaffold is W-linked.
        else if (male_mean <= 0.10 && female_mean >= 0.20) {
            classification = "possible_W_linked"
        }

        # Write one result per scaffold.
        printf "%s\t%s\t%.4f\t%.4f\t%s\t%s\n", \
               scaffold, scaffold_length[scaffold], \
               male_mean, female_mean, ratio, classification
    }
}

# First input: provisional sex assignments.
# Second input: combined normalised scaffold depths.
' "${SEX_RESULTS}" "${ALL_DEPTHS}" |

# Preserve the header and sort the remaining records:
#
#   primary key:   classification, alphabetically
#   secondary key: scaffold length, largest first
{
    IFS= read -r header
    printf '%s\n' "$header"
    sort -t $'\t' -k6,6 -k2,2nr
} > "${SEX_LINKED_RESULTS}"

# Print the final result location to the Slurm log.
echo "Possible sex-linked scaffolds: ${SEX_LINKED_RESULTS}"






exit 0

# Later in R

setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Sex_identification")

library(tidyverse)
library(ggplot2)
library(ggrepel) # For labels
library(readr)

# ============================================================
# INPUT FILES
# ============================================================

sex_file <- "sex_results.tsv"
depth_file <- "all_normalised_depths.tsv"
candidate_file <- "possible_sex_linked_scaffolds.tsv"

# Candidate scaffolds to examine
candidate_scaffolds <- c(
  "scaffold_5",
  "scaffold_53",
  "scaffold_196"
)

# ============================================================
# READ DATA
# ============================================================

sex_results <- read_tsv(
  sex_file,
  show_col_types = FALSE
)

all_depths <- read_tsv(
  depth_file,
  show_col_types = FALSE
)

candidate_results <- read_tsv(
  candidate_file,
  show_col_types = FALSE
)

# Add the predicted sex to the individual-level scaffold-depth data
plot_data <- all_depths %>%
  left_join(
    sex_results %>% select(Sample, Predicted_sex),
    by = "Sample"
  )

# Set a consistent order for sex categories
plot_data <- plot_data %>%
  mutate(
    Predicted_sex = factor(
      Predicted_sex,
      levels = c("male_ZZ", "female_ZW", "ambiguous")
    )
  )




# PLot individual depths
candidate_depths <- plot_data %>%
  filter(
    Scaffold %in% candidate_scaffolds,
    Predicted_sex != "ambiguous"
  ) %>%
  mutate(
    Scaffold = factor(
      Scaffold,
      levels = candidate_scaffolds
    )
  )

candidate_plot <- ggplot(
  candidate_depths,
  aes(
    x = Predicted_sex,
    y = Normalised_depth
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted"
  ) +
  facet_wrap(
    ~ Scaffold,
    nrow = 1
  ) +
  labs(
    x = "Predicted sex",
    y = "Normalised scaffold depth",
    title = "Individual coverage patterns for candidate Z-linked scaffolds",
    subtitle = "Expected Z-linked depth: males ≈ 1 and females ≈ 0.5"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

candidate_plot

ggsave(
  "candidate_Z_scaffold_individual_depths.png",
  candidate_plot,
  width = 11,
  height = 5
)




# HEATMAP CREATION
# ============================================================
# HEATMAP OF NORMALISED DEPTH ACROSS SEX-LINKED AND CONTROL SCAFFOLDS
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

sex_file <- "sex_results.tsv"
depth_file <- "all_normalised_depths.tsv"
candidate_file <- "possible_sex_linked_scaffolds.tsv"


# ------------------------------------------------------------
# Select scaffolds
# ------------------------------------------------------------

selected_scaffolds <- c(
  "scaffold_1",
  "scaffold_5",
  "scaffold_53",
  "scaffold_196"
)

selected_scaffolds

# ------------------------------------------------------------
# Prepare heatmap data
# ------------------------------------------------------------

heatmap_data <- all_depths %>%
  filter(Scaffold %in% selected_scaffolds) %>%
  left_join(
    sex_results %>%
      select(
        Sample,
        Predicted_sex,
        Sex_scaffold_depth = Normalised_depth
      ),
    by = "Sample"
  ) %>%
  mutate(
    Predicted_sex = factor(
      Predicted_sex,
      levels = c("male_ZZ", "female_ZW", "ambiguous")
    )
  )


# ------------------------------------------------------------
# Order samples
# ------------------------------------------------------------

sample_order <- sex_results %>%
  pull(Sample)

heatmap_data <- heatmap_data %>%
  mutate(
    Sample = factor(
      Sample,
      levels = rev(sample_order)
    )
  )


# ------------------------------------------------------------
# Order scaffolds
# ------------------------------------------------------------

heatmap_data <- heatmap_data %>%
  mutate(
    Scaffold = factor(
      Scaffold,
      levels = c(
        "scaffold_1",
        "scaffold_5",
        "scaffold_53",
        "scaffold_196"
      )
    )
  )


# ------------------------------------------------------------
# Plot heatmap
# ------------------------------------------------------------

depth_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Scaffold,
    y = Sample,
    fill = Normalised_depth
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.2
  ) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1.5),
    oob = scales::squish,
    breaks = c(0, 0.5, 1.0, 1.5),
    name = "Normalised\ndepth"
  ) +
  facet_grid(
    Predicted_sex ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  labs(
    x = "Scaffold",
    y = "Sample",
    title = "Normalised scaffold depth across Orange-bellied Parrot samples",
    subtitle = paste(
      "Control autosomal scaffold (scaffold 1) and candidate sex-linked scaffolds;"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.placement = "outside",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

depth_heatmap


ggsave(
  filename = "normalised_depth_sex_scaffold_heatmap.png",
  plot = depth_heatmap,
  width = 11,
  height = 10,
  dpi = 300
)

