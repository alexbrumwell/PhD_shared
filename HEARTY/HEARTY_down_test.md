# HEARTY - Test downsampling 09-06-26

We wanted to test the effects of downsampling using two methods:

1. Downsample a BAM and then run HEARTY
2. Using a downsmpling parameter within HEARTY  (i.e. using the ANGSD min and max depth parameter)

## Workflow

1. Run HEARTY ↓
1. Plot minor allele balance / frequency outputs
   ↓
1. Parse HEARTY het window files and calculate pooled thresholds
   ↓
1. Generate run_roh_sensitivity.sh
   ↓
1. Submit/run run_roh_sensitivity.sh
   ↓
1. Summarise ROH sensitivity outputs into ROH_summary_all.tsv
   ↓
1. Plot ROH/FROH sensitivity from ROH_summary_all.tsv

## Parameters I used for this workflow (ADJUST IF NEEDED)

- minimum depth: 4
- maximum depth: 50
- HEARTY minor allele balance thresholds: 0.15, 0.20, 0.25, 0.30
- region: scaffold_1
- ROH thresholds: 0.0001, 0.0005, 0.001

<br>

## Method A: Normal downsampling BAMs

### 0. Downsample BAMs to 7x

```bash
#!/bin/bash
#SBATCH --job-name=Downsample 
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

timestamp() { date '+%F %T'; }

# === Load modules ===
module load samtools/1.21

# === MODIFY ===
TARGET_COVERAGES=("7" "5")  # target coverages to downsample to

DIR_BAM="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP" # Direcotry where BAM files are located
LIST_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected.list" # List of sample names (without .rmdup.bam extension) to process

OUTDIR="${DIR_BAM}/Downsampled_normal"  # Outdirectory where downsampled BAM files will be saved
OUTDIR_7X="${OUTDIR}/7x"
OUTDIR_5X="${OUTDIR}/5x"
mkdir -p "$OUTDIR" "$OUTDIR_7X" "$OUTDIR_5X"


# ===== Output summary files =====
LIST_COVERAGE="${OUTDIR}/coverage_original.tsv" # Stores original coverage for each input BAM

QC_OUT="${OUTDIR}/downsampling_qc.tsv" # Stores downsampling information and observed coverage after downsampling

LIST_DOWN7X="${OUTDIR}/Down_7x.list" # Lists BAMs to use for the 7x and 5x datasets
LIST_DOWN5X="${OUTDIR}/Down_5x.list"

# Empty old output files before starting
: > "$LIST_COVERAGE"
: > "$QC_OUT"
: > "$LIST_DOWN7X"
: > "$LIST_DOWN5X"

# Add header to the QC table
echo -e "Sample\tInput_BAM\tOriginal_coverage\tTarget_coverage\tFraction\tOutput_BAM\tObserved_coverage" > "$QC_OUT"


# ===== Function to calculate mean genome-wide depth =====
calc_cov() {
    local bam="$1"

    # samtools depth -a reports depth at every position, including zero-depth sites.
    # awk sums depth across all positions and divides by the number of positions.
    samtools depth -@ "$SLURM_CPUS_PER_TASK" -a "$bam" | awk '{sum += $3; n++} END {if (n > 0) print sum/n; else print "NA"}'
}

# ===== Function to downsample one BAM to one target coverage =====
downsample_if_needed() {
    local bam="$1"        # Input BAM
    local sample="$2"     # Sample name
    local coverage="$3"   # Original observed coverage
    local target="$4"     # Target coverage, e.g. 7 or 5
    local outdir="$5"     # Output directory
    local list_out="$6"   # File to store output BAM path

    local output_bam="${outdir}/${sample}.down${target}x.bam"

    # Only downsample if original coverage is greater than the target coverage
    if awk -v c="$coverage" -v t="$target" 'BEGIN {exit !(c > t)}'; then

        # Calculate the fraction of reads to keep
        # Example: if coverage is 10x and target is 5x, fraction = 5/10 = 0.5
        frac=$(awk -v t="$target" -v c="$coverage" 'BEGIN {print t/c}')

        echo "[$(timestamp)] Downsampling ${sample}: ${coverage}x -> ${target}x; fraction=${frac}"

        # Only run downsampling if output BAM does not already exist
        if [[ ! -f "$output_bam" ]]; then

            # Downsample BAM using samtools view.
            # -b writes BAM output.
            # -s keeps a random fraction of reads.
            samtools view \
                -@ "$SLURM_CPUS_PER_TASK" \
                -b \
                -s "${frac}" \
                "$bam" \
                -o "$output_bam"

            # Index the downsampled BAM
            samtools index -@ "$SLURM_CPUS_PER_TASK" "$output_bam"
						samtools stats -@ "$SLURM_CPUS_PER_TASK" $output_bam > ${output_bam%.bam}.stats.txt
        else
            echo "[$(timestamp)] Output already exists: $output_bam"
        fi

        # Calculate observed coverage after downsampling
        observed=$(calc_cov "$output_bam")

        # Add downsampled BAM path to the relevant BAM list
        echo "$output_bam" >> "$list_out"

    else
        # If original coverage is already below the target, keep the original BAM
        echo "[$(timestamp)] ${sample} has coverage ${coverage}x, which is <= ${target}x. Keeping original BAM."

        output_bam="$bam"
        frac="NA"
        observed="$coverage"

        # Add original BAM path to the relevant BAM list
        echo "$output_bam" >> "$list_out"
    fi

    # Write one line to the QC summary table
    echo -e "${sample}\t${bam}\t${coverage}\t${target}\t${frac}\t${output_bam}\t${observed}" >> "$QC_OUT"
}


# ===== Main loop over samples =====
while read -r sample; do

    # Skip empty lines in the sample list
    [[ -z "$sample" ]] && continue

    # Construct the expected BAM path
    bam="${DIR_BAM}/${sample}.rmdup.bam"

    # Check that the BAM exists before processing
    if [[ ! -f "$bam" ]]; then
        echo "Error: BAM file not found: $bam" >&2
        continue
    fi

    # Calculate original coverage
    echo "[$(timestamp)] Calculating coverage for $sample..."
    coverage=$(calc_cov "$bam")

    # Save original coverage
    echo -e "${sample}\t${bam}\t${coverage}" >> "$LIST_COVERAGE"

    # Downsample to 7x and 5x
    downsample_if_needed "$bam" "$sample" "$coverage" "7" "$OUTDIR_7X" "$LIST_DOWN7X"
    downsample_if_needed "$bam" "$sample" "$coverage" "5" "$OUTDIR_5X" "$LIST_DOWN5X"

done < "$LIST_SAMPLES"

# ===== Final log messages =====
echo "[$(timestamp)] Downsampling complete."
echo "Original coverage table: $LIST_COVERAGE"
echo "Downsampling QC table: $QC_OUT"
echo "7x BAM list: $LIST_DOWN7X"
echo "5x BAM list: $LIST_DOWN5X"
```

### 1. Running HEARTY python script

First we need to actually run the HEARTY python script with the parameters we have set earlier.

```bash
#!/bin/bash
#SBATCH --job-name=HEARTY_down7x
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-24%20
#SBATCH --output=jobs/%x_%A_%a.out
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk
 
# Script to run HEARTY to test allele balances per sample

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21 python

DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/HEARTY/down_normal"
cd $DIR

BAM_LIST="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled_normal/Down_7x.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
#LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"

# Get sample name from BAM list based on SLURM_ARRAY_TASK_ID
BAM=$(cat $BAM_LIST | grep "Downsampled" | sed -n "${SLURM_ARRAY_TASK_ID}p" )
if [ -z "${BAM}" ]; then
    echo "No BAM file found for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

SAMPLE=$(basename "${BAM}" .down7x.bam)

echo "Processing sample ${SAMPLE} with BAM file ${BAM}"

# Apply HEARTY

#All scaffold 1
MAXDEPTH="50"
MINDEPTH="4"

python3 /projects/erode/people/qhk515/software/HEARTY/HEARTY.py --bam ${BAM} \
    -o ${SAMPLE}.allscaf1.downnormal --cores ${SLURM_CPUS_PER_TASK} \
    --regions "scaffold_1" --angsd /opt/shared_software/shared_envmodules/conda/angsd-0.940/bin/angsd --force \
    --threshold 0.15 \
    --threshold 0.20 \
    --threshold 0.25 \
    --threshold 0.30 \
    --min-depth ${MINDEPTH} \
    --max-depth ${MAXDEPTH} \
    --reuse-existing
echo "Finished running HEARTY for $SAMPLE "

exit 0
```

> Continue steps in the normal HEARTY pipeline (see other doc) but using these output files.
>
>
>
>



#### R script to plot

```bash
# Script to plot frequency output from HEARTY

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)


# OPTIONAL:
# To join pdfs together
library(pdftools)
library(stringr)

pdfs <- list.files(
  #path = "path/to/pdfs",
  pattern = "minorfreq_overlay\\.pdf$",
  full.names = TRUE
)

# Natural sorting (important for sample IDs)
pdfs <- str_sort(pdfs, numeric = TRUE)

pdf_combine(pdfs, output = "combined_minorfreq_plots.pdf")



# ------------------------------------------------------------
# 1. Directory containing HEARTY minorfreq files
# ------------------------------------------------------------

setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Het/HEARTY/down_cap/minorfreq")

plot_dir <- "plots"

dir.create(plot_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Find raw minorfreq files
#    Adjust pattern if needed
# ------------------------------------------------------------

files <- list.files(
  pattern = "minorfreq\\.txt$",
  full.names = TRUE
)

# Exclude summary files
files <- files[!str_detect(files, "summary")]
#files <- files[str_detect(files, "downcap")]

files

# ------------------------------------------------------------
# 3. Read and combine all individuals
# ------------------------------------------------------------

pair_types <- tibble(
  pair = c("AC", "AG", "AT", "CG", "CT", "GT"),
  mutation_type = c(
    "Transversion",
    "Transition",
    "Transversion",
    "Transversion",
    "Transition",
    "Transversion"
  )
)

minor_df <- map_dfr(files, function(f) {
  
  sample_name <- basename(f) %>%
    str_remove("\\..*$")
  
  read_tsv(f, show_col_types = FALSE) %>%
    pivot_longer(
      cols = c(AC, AG, AT, CG, CT, GT),
      names_to = "pair",
      values_to = "count"
    ) %>%
    left_join(pair_types, by = "pair") %>%
    mutate(
      sample = sample_name,
      file = basename(f)
    )
})

minor_df <- minor_df %>%
  mutate(
    pair = factor(pair, levels = c("AC", "AG", "AT", "CG", "CT", "GT")),
    mutation_type = factor(mutation_type, levels = c("Transition", "Transversion"))
  ) %>%
  mutate(
    Group = if_else(
      str_starts(sample, "O-"),
      "Modern",
      "Historical"
    )
  )

str(minor_df)

# Samples used for downsampling
down_samples <- c(
  "O-76339-001",
  "O-76349-001",
  "O-76521-001",
  "O-76523-001",
  "O-76526-001",
  "O-76538-001",
  "O-81171-001",
  "O-81175-001",
  "O-81179-001",
  "O-81182-001",
  "O-81186-001",
  "O-81188-001",
  "O-81192-001",
  "O-81193-001",
  "O-81196-001",
  "O-81201-001",
  "O-81203-001",
  "OBPI_501",
  "OBPI_502",
  "OBPI_503",
  "OBPI_504",
  "OBPII_ReHoI_504",
  "OBPII_ReHoI_507",
  "OBPII_ReHoI_508"
)

minor_df <- minor_df %>%
  filter(sample %in% down_samples)

# Add coverage
coverage_df <- tribble(
  ~sample, ~MeanCoverage,
  "O-76339-001", 7.18,
  "O-76349-001", 7.43,
  "O-76350-001", 6.95,
  "O-76521-001", 7.04,
  "O-76523-001", 7.08,
  "O-76526-001", 7.07,
  "O-76538-001", 7.07,
  "O-81171-001", 7.10,
  "O-81175-001", 7.09,
  "O-81179-001", 7.26,
  "O-81182-001", 6.97,
  "O-81186-001", 6.84,
  "O-81188-001", 7.01,
  "O-81192-001", 7.16,
  "O-81193-001", 7.10,
  "O-81196-001", 6.95,
  "O-81198-001", 6.08,
  "O-81201-001", 7.01,
  "O-81203-001", 6.93,
  "OBPI_501", 24.06,
  "OBPI_502", 9.97,
  "OBPI_503", 7.24,
  "OBPI_504", 17.67,
  "OBPI_505", 2.42,
  "OBPI_506", 4.86,
  "OBPI_507", 5.77,
  "OBPI_508", 3.91,
  "OBPII_ReHoI_501", 5.10,
  "OBPII_ReHoI_502", 6.71,
  "OBPII_ReHoI_503", 3.94,
  "OBPII_ReHoI_504", 7.60,
  "OBPII_ReHoI_505", 2.64,
  "OBPII_ReHoI_506", 4.03,
  "OBPII_ReHoI_507", 9.06,
  "OBPII_ReHoI_508", 10.24
)

# Join coverage
minor_df <- minor_df %>%
  left_join(coverage_df, by = "sample")


#Coverage labels for plots
minor_df <- minor_df %>%
  mutate(sample_cov = paste0(sample, " (", MeanCoverage, "x)"))
# Reorder samples
minor_df <- minor_df %>%
  mutate(sample_cov = forcats::fct_reorder(sample_cov, MeanCoverage))

# ------------------------------------------------------------
# 4. Plot: one panel per individual, shared y axis
# ------------------------------------------------------------

# PLOT 1: raw frequencies
p <- ggplot(
  data = minor_df %>% filter(Group == "Historical" | sample ==  "O-76339-001") %>% filter(sample != "OBPI_504"),
  #data = minor_df %>% filter(sample ==  "O-76339-001"),
  aes(
    x = freq,
    y = count,
    colour = pair,
    shape = mutation_type
  )
) +
  geom_point(size = 0.5, alpha = 0.2) +
  geom_line(aes(group = pair), linewidth = 0.5, alpha = 0.7) +
  facet_wrap(Group ~ sample_cov, scales = "free", ncol = 4) +
  theme_bw() +
  labs(
    x = "Allele balance / minor allele frequency",
    y = "Number of sites",
    colour = "Base change",
    shape = "Mutation type",
    title = "Minor allele balance distribution across historical samples OBP individuals",
    subtitle = "METHOD B: HEARTY DOWNSAMPLING. Samples ordered by increasing coverage"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
p
ggsave(
  filename = file.path(plot_dir, "01_OBP_minorfreq_all_individuals_diff_y.pdf"),
  plot = p,
  width = 12,
  height = 8
)

 
# PLOT 2: Cumulative thresholds
# ------------------------------------------------------------
# Cumulative threshold plot: different mutation pairs
# ------------------------------------------------------------

cdf_pair_df <- minor_df %>%
  group_by(sample, Group, sample_cov, pair, mutation_type, freq) %>%
  summarise(
    total_sites = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample, Group, sample_cov, pair, mutation_type) %>%
  arrange(freq, .by_group = TRUE) %>%
  mutate(
    cum_sites = cumsum(total_sites),
    cum_prop_sites = cum_sites / sum(total_sites)
  ) %>%
  ungroup()

p_cdf_pair <- ggplot(
  data = cdf_pair_df %>%
    filter(Group == "Historical" | sample == "O-76339-001") %>%
    filter(sample != "OBPI_504"),
  aes(x = freq, y = cum_prop_sites, colour = pair, linetype = mutation_type)
) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(
    xintercept = c(0.05, 0.10, 0.20, 0.25, 0.30),
    linetype = "dashed",
    alpha = 0.35
  ) +
  facet_wrap(~ sample_cov, scales = "fixed", ncol = 4) +
  theme_bw() +
  labs(
    x = "Allele balance / minor allele frequency threshold",
    y = "Cumulative proportion of sites",
    colour = "Base change",
    linetype = "Mutation type",
    title = "Cumulative allele balance distribution by mutation class",
    subtitle = "Each line shows the cumulative proportion of sites for one base-change class"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_cdf_pair

ggsave(
  filename = file.path(plot_dir, "04_OBP_cumulative_thresholds_by_mutation_pair.pdf"),
  plot = p_cdf_pair,
  width = 12,
  height = 8
)


# PLOT: TRANSITION/ TRANSVERSION RATIO
titv_df <- minor_df %>%
  
  group_by(Group, sample, freq, mutation_type) %>%
  
  summarise(count = sum(count), .groups = "drop") %>%
  
  pivot_wider(
    names_from = mutation_type,
    values_from = count
  ) %>%
  
  mutate(
    TiTv = Transition / Transversion
  )

p3 <- ggplot(
  titv_df,
  aes(x = freq, y = TiTv, group = sample, colour = Group)
) +
  geom_line(alpha = 0.25) +
  stat_summary(
    aes(group = Group),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  theme_bw() +
  labs(
    x = "Allele balance",
    y = "Transition / Transversion ratio",
    title = "Ti/Tv ratio across allele balance"
  )
p3 
ggsave(file.path(plot_dir, "03_TSTV_hist_mod.png"), p3, width = 10, height = 8)


```





## Method B: Downsample using HEARTYdepth caps

### 1. Run HEARTY script

This is basically the same script as the original but an added parameter of `--downsample-depth 7`

```bash
#!/bin/bash
#SBATCH --job-name=HEARTY_down7x_cap
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%20
#SBATCH --output=jobs/%x_%A_%a.out
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk
 
# Script to run HEARTY to test allele balances per sample

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21 python

DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/HEARTY/down_cap"
cd $DIR

BAM_LIST="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_OBPs_combined.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
#LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"

# Get sample name from BAM list based on SLURM_ARRAY_TASK_ID
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BAM_LIST})
if [ -z "${BAM}" ]; then
    echo "No BAM file found for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


SAMPLE=$(basename "${BAM}" .rmdup.bam)

echo "Processing sample ${SAMPLE} with BAM file ${BAM}"

# Apply HEARTY

#All scaffold 1
MAXDEPTH="50"
MINDEPTH="4"

python3 /projects/erode/people/qhk515/software/HEARTY/HEARTY.py --bam ${BAM} \
    -o ${SAMPLE}.allscaf1.downcap --cores ${SLURM_CPUS_PER_TASK} \
    --regions "scaffold_1" --angsd /opt/shared_software/shared_envmodules/conda/angsd-0.940/bin/angsd --force \
    --threshold 0.15 \
    --threshold 0.20 \
    --threshold 0.25 \
    --threshold 0.30 \
    --min-depth ${MINDEPTH} \
    --max-depth ${MAXDEPTH} \
		--downsample-depth 7 --downsample-reps 10 

echo "Finished running HEARTY for $SAMPLE "

```

> Continue steps in the normal HEARTY pipeline (see other doc) but using these output files.