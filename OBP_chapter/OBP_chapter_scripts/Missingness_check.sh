#!/bin/bash
#SBATCH --job-name=Missingness_check
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
 #SBATCH --partition=cpuqueue
 #SBATCH --qos=normal
 #SBATCH --output=jobs/%x_%j.log
 #SBATCH --error=jobs/%x_%j.err
 #SBATCH --mail-type=END,FAIL
 #SBATCH --mail-user=amb230@kent.ac.uk

 set -euo pipefail

 module purge
 module load perl/5.38.0 gsl/2.5 bcftools/1.21 gcc/13.2.0 openjdk/20.0.0 vcftools R

 timestamp() { date '+%F %T'; }

 # === Directories ===
 INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_GroupMissFiltered"
 OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_GroupMissFiltered/missingness_check"
 PLOTDIR="${OUTDIR}/plots"
 mkdir -p "${OUTDIR}" "${PLOTDIR}"

 CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

 echo "$(timestamp) Starting per-sample missingness analysis"

 # === Loop through all chromosomes ===
 while IFS= read -r chrom; do
 INVCF="${INDIR}/OBP_merged.26.${chrom}.groupmiss80.vcf.gz"

 if [[ ! -s "${INVCF}" ]]; then
 echo "WARNING: Input VCF not found: ${INVCF}"
 continue
 fi

 echo "$(timestamp) Calculating per-sample missingness for ${chrom}"

 # Calculate per-individual missingness
 vcftools \
 --gzvcf "${INVCF}" \
 --missing-indv \
 --out "${OUTDIR}/OBP.${chrom}.groupmiss80"

 echo "$(timestamp) Completed for ${chrom}"
 done < "${CHROM_LIST}"

 echo "$(timestamp) All chromosomes processed. Generating plots..."

 # === Generate plots with R ===
 Rscript --vanilla - "${OUTDIR}" "${PLOTDIR}" <<'RSCRIPT'
 library(data.table)
 library(ggplot2)
 library(dplyr)

 args <- commandArgs(trailingOnly = TRUE)
 input_dir <- args[1]
 output_dir <- args[2]

 if (!dir.exists(output_dir)) {
 dir.create(output_dir, recursive = TRUE)
 }

 # Find all .imiss files
 imiss_files <- list.files(input_dir, pattern = "\\.imiss$", full.names = TRUE)

 if (length(imiss_files) == 0) {
 stop("No .imiss files found in ", input_dir)
 }

 cat("Found", length(imiss_files), "per-chromosome missingness files\n")

 # Read and combine all files
 all_data <- data.table()

 for (file in imiss_files) {
 chrom <- gsub(".*OBP\\.(.*)\\.groupmiss80\\.imiss", "\\1", file)
 data <- fread(file)
 data$CHROM <- chrom
 all_data <- rbind(all_data, data)
 }

 # Calculate per-sample statistics across all chromosomes
sample_stats <- all_data %>%
group_by(INDV) %>%
summarise(
 mean_F_MISS = mean(F_MISS, na.rm = TRUE),
 sd_F_MISS = sd(F_MISS, na.rm = TRUE),
 min_F_MISS = min(F_MISS, na.rm = TRUE),
 max_F_MISS = max(F_MISS, na.rm = TRUE),
 n_chroms = n(),
 n_missing_chroms = sum(is.na(F_MISS)),   # worth tracking
 .groups = 'drop'
 ) %>%
 arrange(desc(mean_F_MISS))

 # Identify outliers (mean + 2*SD)
 overall_mean <- mean(sample_stats$mean_F_MISS)
 overall_sd <- sd(sample_stats$mean_F_MISS)
 threshold <- overall_mean + 2 * overall_sd

 sample_stats$is_outlier <- sample_stats$mean_F_MISS > threshold

 # Write summary table
 write.csv(sample_stats,
 file.path(output_dir, "sample_missingness_summary.csv"),
 row.names = FALSE)

 # Write full data table (all samples x all chromosomes)
 write.csv(all_data,
 file.path(output_dir, "sample_missingness_full_data.csv"),
 row.names = FALSE)

 cat("\nPer-sample missingness summary:\n")
 print(sample_stats)

 cat("\n\nOutlier samples (mean  mean + 2*SD):\n")
 outliers <- sample_stats %>% filter(is_outlier)
 if (nrow(outliers) > 0) {
 print(outliers)
 } else {
 cat("No outliers detected\n")
 }

 # === Plot 1: Bar plot of mean missingness per sample ===
 p1 <- ggplot(sample_stats %>% arrange(mean_F_MISS),
 aes(x = reorder(INDV, mean_F_MISS), y = mean_F_MISS, fill = is_outlier)) +
 geom_bar(stat = "identity") +
 geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 1) +
 geom_hline(yintercept = overall_mean, linetype = "dashed", color = "blue", size = 1) +
 scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "orangered")) +
 theme_minimal() +
 theme(
 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
 panel.grid.major.y = element_line(color = "grey90"),
 legend.position = "top"
 ) +
 labs(
 title = "Per-Sample Missingness (80% filter)",
 x = "Sample ID",
 y = "Mean Proportion Missing",
 fill = "Outlier",
 subtitle = paste0("Blue line = mean (", round(overall_mean, 4),
 "), Red line = threshold (", round(threshold, 4), ")")
 )

 ggsave(file.path(output_dir, "sample_missingness_barplot.png"),
 p1, width = 14, height = 6, dpi = 300)

 # === Plot 2: Box plot showing per-sample variation across chromosomes ===
 all_data_melted <- all_data %>%
 mutate(INDV = factor(INDV)) %>%
 arrange(INDV)

# Calculate mean per sample for ordering
 sample_order <- all_data_melted %>%
 group_by(INDV) %>%
 summarise(mean_miss = mean(F_MISS)) %>%
 arrange(mean_miss) %>%
 pull(INDV)

 all_data_melted$INDV <- factor(all_data_melted$INDV, levels = sample_order)

> p2 <- ggplot(all_data_melted, aes(x = INDV, y = F_MISS, fill = INDV)) +
 geom_boxplot(show.legend = FALSE) +
 theme_minimal() +
 theme(
 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
 panel.grid.major.y = element_line(color = "grey90")
 ) +
 labs(
 title = "Distribution of Missingness Across Chromosomes",
 x = "Sample ID",
 y = "Proportion Missing"
 )

 ggsave(file.path(output_dir, "sample_missingness_boxplot.png"),
 p2, width = 14, height = 6, dpi = 300)

 # === Plot 3: Scatter plot of mean vs SD ===
 p3 <- ggplot(sample_stats, aes(x = mean_F_MISS, y = sd_F_MISS, color = is_outlier, label = INDV)) +
 geom_point(size = 3, alpha = 0.7) +
 scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "orangered")) +
 theme_minimal() +
 theme(
 panel.grid.major = element_line(color = "grey90"),
 legend.position = "top"
 ) +
 labs(
 title = "Mean vs. SD of Missingness per Sample",
 x = "Mean Proportion Missing",
 y = "SD of Proportion Missing",
 color = "Outlier"
 )

 # Add labels to outliers
> if (nrow(outliers) > 0) {
> p3 <- p3 + geom_text(data = outliers, aes(label = INDV),
> hjust = -0.2, vjust = -0.2, size = 3)
> }
>
 ggsave(file.path(output_dir, "sample_missingness_scatter.png"),
 p3, width = 10, height = 8, dpi = 300)

 cat("\n\nPlots saved to:", output_dir, "\n")
 cat("Files created:\n")
 cat("  - sample_missingness_summary.csv\n")
 cat("  - sample_missingness_barplot.png\n")
 cat("  - sample_missingness_boxplot.png\n")
 cat("  - sample_missingness_scatter.png\n")
 RSCRIPT

 echo "$(timestamp) Analysis complete!"