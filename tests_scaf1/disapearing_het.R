# ============================================================
# Compare ANGSD GL het sites against SNPAD, BCFTOOLS and ANGSD_VCF
# ============================================================

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(writexl)
library(patchwork)
library(vroom)

# Optional
# library(httpgd)
# hgd()

setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Het/tests_scaf1_angsdvssnpad/disappearing_het_sites")

outdir <- "plots_angsd_GL_vs_all_callers"
dir.create(outdir, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Read only selected sample files
# ------------------------------------------------------------

samples_to_keep <- c(
  "O-76339-001",
  "OBPI_501","OBPI_502","OBPI_503","OBPI_504",
  "OBPI_505","OBPI_506","OBPI_507","OBPI_508",
  "OBPII_ReHoI_501","OBPII_ReHoI_502","OBPII_ReHoI_503","OBPII_ReHoI_504",
  "OBPII_ReHoI_505","OBPII_ReHoI_506","OBPII_ReHoI_507","OBPII_ReHoI_508"
  )
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

files_to_read <- paste0(
  samples_to_keep,
  "_final_angsd_GL_hets_vs_all_vcf_callers.tsv"
)

# Check files exist
missing_files <- files_to_read[!file.exists(files_to_read)]
if (length(missing_files) > 0) {
  stop("These files are missing:\n", paste(missing_files, collapse = "\n"))
}

df <- vroom(
  files_to_read,
  col_types = cols(.default = "c"),
  #id = "sample"
)
names(df)


# ------------------------------------------------------------
# 2. Clean duplicated column names
# ------------------------------------------------------------

# Because the final file joins ANGSD evidence + VCF result,
# readr may create names like key...1, chrom...3, key...21 etc.
# This tries to standardise the important columns.
names(df)[1:37] <- c(
  "key_angsd",
  "sample",
  "chrom_angsd",
  "pos_angsd",
  "ANGSD_A1",
  "ANGSD_A2",
  "ANGSD_pAA",
  "ANGSD_pAB",
  "ANGSD_pBB",
  "ANGSD_call",
  "maf_ref",
  "maf_major",
  "maf_minor",
  "ANGSD_count_chr",
  "ANGSD_count_pos",
  "ANGSD_totDepth",
  "ANGSD_totA",
  "ANGSD_totC",
  "ANGSD_totG",
  "ANGSD_totT",
  "vcf_key",
  "vcf_chrom",
  "vcf_pos",
  "vcf_REF",
  "vcf_ALT",
  "method",
  "present",
  "GT",
  "DP",
  "allele_A_raw",
  "allele_C_raw",
  "allele_G_raw",
  "allele_T_raw",
  "REF_reads",
  "ALT_reads",
  "GQ_or_GP",
  "outcome"
)


# ------------------------------------------------------------
# 4. Clean and derive useful columns
# ------------------------------------------------------------
sum_comma_counts <- function(x) {
  sapply(strsplit(as.character(x), ","), function(y) {
    y <- suppressWarnings(as.numeric(y))
    if (all(is.na(y))) NA_real_ else sum(y, na.rm = TRUE)
  })
}
allele_map <- c("0" = "A", "1" = "C", "2" = "G", "3" = "T")

df$GQ_or_GP <- as.numeric(df$GQ_or_GP)

# Determine outcome
classify_vcf_outcome <- function(present, gt) {
  present <- as.character(present)
  gt <- as.character(gt)
  
  gt <- stringr::str_replace(gt, ":.*$", "")
  
  a1 <- stringr::str_split_fixed(gt, "[/|]", 2)[, 1]
  a2 <- stringr::str_split_fixed(gt, "[/|]", 2)[, 2]
  
  dplyr::case_when(
    present == "no" ~ "absent_from_vcf",
    
    is.na(gt) | gt %in% c("NA", ".", "./.", ".|.") ~ "missing_genotype",
    
    a1 == "." | a2 == "." ~ "missing_genotype",
    
    a1 == "0" & a2 == "0" ~ "converted_to_hom_ref",
    
    a1 == a2 & a1 != "0" ~ "converted_to_hom_alt",
    
    a1 != a2 ~ "retained_het",
    
    TRUE ~ "other"
  )
}


# 1. Basic types only
df2 <- df %>%
  mutate(
    sample = as.character(sample),
    method = factor(method, levels = c("SNPAD", "BCFTOOLS", "ANGSD_VCF"))
  )

# Join coverage
df2 <- df2 %>%
  left_join(coverage_df, by = "sample")

stopifnot(!any(is.na(df2$MeanCoverage)))

# 2. Decode ANGSD alleles
df2 <- df2 %>%
  mutate(
    ANGSD_A1_base = allele_map[as.character(ANGSD_A1)],
    ANGSD_A2_base = allele_map[as.character(ANGSD_A2)]
  )


# 3. Convert numeric columns
numeric_cols <- c(
  "pos_angsd", "vcf_pos",
  "ANGSD_A1", "ANGSD_A2",
  "ANGSD_pAA", "ANGSD_pAB", "ANGSD_pBB",
  "ANGSD_totDepth", "ANGSD_totA", "ANGSD_totC",
  "ANGSD_totG", "ANGSD_totT",
  "DP"
)

df2 <- df2 %>%
  mutate(across(all_of(numeric_cols), as.numeric))

# 4. Parse VCF allele counts
df2 <- df2 %>%
  mutate(
    allele_A = sum_comma_counts(allele_A_raw),
    allele_C = sum_comma_counts(allele_C_raw),
    allele_G = sum_comma_counts(allele_G_raw),
    allele_T = sum_comma_counts(allele_T_raw)
  )
# 5. Clean GT
df2 <- df2 %>%
  mutate(
    GT_clean = sub(":.*$", "", as.character(GT))
  )
# 6. Classify outcome
df2 <- df2 %>%
  mutate(
    outcome_clean = classify_vcf_outcome(
      present = present,
      gt = GT_clean
    )
  )
# 7. Apply factor levels
df2 <- df2 %>%
  mutate(
    outcome_clean = factor(
      outcome_clean,
      levels = c(
        "retained_het",
        "converted_to_hom_ref",
        "converted_to_hom_alt",
        "missing_genotype",
        "absent_from_vcf",
        "other"
      )
    )
  )

stopifnot(
  "present" %in% names(df),
  "GT" %in% names(df),
  "method" %in% names(df)
)

df %>%
  count(method, present, GT, sort = TRUE) %>%
  print(n = 50)

# ------------------------------------------------------------
# 5. ANGSD read-balance metrics
# ------------------------------------------------------------
df2 <- df2 %>%
  mutate(
    ANGSD_A1_reads = case_when(
      ANGSD_A1_base == "A" ~ ANGSD_totA,
      ANGSD_A1_base == "C" ~ ANGSD_totC,
      ANGSD_A1_base == "G" ~ ANGSD_totG,
      ANGSD_A1_base == "T" ~ ANGSD_totT,
      TRUE ~ NA_real_
    ),

    ANGSD_A2_reads = case_when(
      ANGSD_A2_base == "A" ~ ANGSD_totA,
      ANGSD_A2_base == "C" ~ ANGSD_totC,
      ANGSD_A2_base == "G" ~ ANGSD_totG,
      ANGSD_A2_base == "T" ~ ANGSD_totT,
      TRUE ~ NA_real_
    ),

    ANGSD_minor_reads = pmin(ANGSD_A1_reads, ANGSD_A2_reads, na.rm = FALSE),
    ANGSD_major_reads = pmax(ANGSD_A1_reads, ANGSD_A2_reads, na.rm = FALSE),
    ANGSD_allelic_depth = ANGSD_A1_reads + ANGSD_A2_reads,

    ANGSD_minor_balance = if_else(
      ANGSD_allelic_depth > 0,
      ANGSD_minor_reads / ANGSD_allelic_depth,
      NA_real_
    ),

    high_conf_angsd_het = ANGSD_pAB >= 0.95,
    balanced_angsd_het = ANGSD_minor_balance >= 0.30,
    high_conf_balanced_het = high_conf_angsd_het & balanced_angsd_het
  )

# ------------------------------------------------------------
# 6. Main summaries
# ------------------------------------------------------------

summary_outcomes <- df2 %>%
  count(sample, MeanCoverage, method, outcome_clean, name = "n_sites") %>%
  group_by(sample, MeanCoverage, method) %>%
  mutate(prop = n_sites / sum(n_sites)) %>%
  ungroup()

summary_by_method <- df2 %>%
  count(sample, MeanCoverage, method, outcome_clean) %>%
  group_by(sample, MeanCoverage, method) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

summary_mechanism <- df2 %>%
  group_by(method, outcome_clean) %>%
  summarise(
    n_sites = n(),
    median_ANGSD_pAB = median(ANGSD_pAB, na.rm = TRUE),
    median_ANGSD_depth = median(ANGSD_totDepth, na.rm = TRUE),
    median_ANGSD_minor_balance = median(ANGSD_minor_balance, na.rm = TRUE),
    median_VCF_DP = median(DP, na.rm = TRUE),
    median_GQ_or_GP = median(GQ_or_GP, na.rm = TRUE),
    prop_high_conf_balanced = mean(high_conf_balanced_het, na.rm = TRUE),
    .groups = "drop"
  )

summary_highconf <- df2 %>%
  filter(high_conf_balanced_het) %>%
  count(sample, method, outcome_clean, name = "n_sites") %>%
  group_by(sample, method) %>%
  mutate(prop = n_sites / sum(n_sites)) %>%
  ungroup()

write_xlsx(
  list(
    outcome_by_sample = summary_outcomes,
    outcome_by_method = summary_by_method,
    mechanism_by_method = summary_mechanism,
    high_conf_balanced = summary_highconf
  ),
  file.path(outdir, "summary_tables_all_callers.xlsx")
)

#Coverage labels for plots
summary_outcomes <- summary_outcomes %>%
  mutate(sample_cov = paste0(sample, " (", MeanCoverage, "x)"))

retention <- summary_outcomes %>%
  filter(outcome_clean == "retained_het") %>%
  mutate(prop_pct = prop * 100)


# PLOT 1: Where do ANGSD GL heterozygous sites go?
# Question: When ANGSD GL calls a site heterozygous, do SNPAD and BCFTOOLS retain it, convert it to hom-ref/hom-alt, or miss it?
summary_outcomes <- summary_outcomes %>%
  mutate(prop_pct = prop * 100)

p1 <- ggplot(summary_outcomes,
             aes(x = forcats::fct_reorder(sample_cov, MeanCoverage),
                 y = prop_pct,
                 fill = outcome_clean)) +
  geom_col() +
  geom_text(
    aes(label = ifelse(prop >= 0.05, scales::percent(prop, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  coord_flip() +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    x = "Sample and mean coverage",
    y = "Percentage of ANGSD GL heterozygous sites",
    fill = "Outcome",
    title = "Outcome of ANGSD GL heterozygous sites per sample",
    subtitle = "Samples ordered by coverage. Percentages per outcome shown if >=5x."
  ) 
print(p1)
ggsave(file.path(outdir, "01_overall_outcome_by_caller.png"),
       p1, width = 14, height = 9, dpi = 300)


# PLOT 2: Which caller retains more ANGSD hets?
# Question: Does SNPAD or BCFTOOLS retain a higher proportion of ANGSD GL heterozygous sites?
p2 <- ggplot(retention,
             aes(x = forcats::fct_reorder(sample_cov, MeanCoverage), y = prop_pct, colour = method, group = method)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Sample mean coverage",
    y = "Percentage retained as heterozygous",
    colour = "Variant caller",
    title = "Which caller retains ANGSD GL heterozygous sites?",
    subtitle = "Samples ordered by coverage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom"
  )

print(p2)
ggsave(file.path(outdir, "02_retained_het_proportion_per_sample.png"),
       p2, width = 8, height = 9, dpi = 300)



#PLOT 3: Are lost ANGSD hets mainly low-depth sites?
# Question: Are ANGSD hets rejected by VCF callers concentrated at low read depth?
p_depth <- ggplot(df2,
                  aes(x = outcome_clean, y = ANGSD_totDepth)) +
  geom_violin(fill = "grey85", colour = NA, scale = "width") +
  stat_summary(fun = median, geom = "point", size = 2, colour = "red") +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    x = "VCF outcome",
    y = "ANGSD total depth",
    title = "Depth distribution of ANGSD GL hets by outcome"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 20))

print(p_depth)

ggsave(file.path(outdir, "03_coverage_het.png"), p_depth, width = 11, height = 6, dpi = 300)


# PLOT 4: Are lost ANGSD hets poorly allele-balanced?
# Question: Do sites rejected by SNPAD/BCFTOOLS have weak minor allele support?
library(ggridges)

p_balance <- ggplot(df2,
                    aes(x = ANGSD_minor_balance,
                        y = outcome_clean,
                        fill = outcome_clean)) +
  geom_density_ridges(scale = 0.9,
                      alpha = 0.7,
                      panel_scaling = FALSE) +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    x = "ANGSD minor allele balance",
    y = "VCF outcome",
    title = "Distribution of allele balance by outcome"
  ) +
  theme(legend.position = "none")
  
print(p_balance)

ggsave(file.path(outdir, "04_balance_het.png"),
       p_balance, width = 11, height = 6, dpi = 300)


# PLOT 5: Does ANGSD give high posterior support to sites rejected by VCF callers?
# Question: Are discordant sites genuinely high-confidence ANGSD hets, or only marginal GL calls?
p_pAB <- ggplot(df2,
                aes(x = outcome_clean, y = ANGSD_pAB)) +
  geom_boxplot(outlier.alpha = 0.1, na.rm = TRUE) +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    x = "VCF outcome",
    y = "ANGSD P(heterozygote)",
    title = "Does ANGSD strongly support hets rejected by VCF callers?"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") 

print(p_pAB)

ggsave(file.path(outdir, "05_posterior_het.png"),
       p_pAB, width = 11, height = 6, dpi = 300)

# Coverage specific plot
p_cov_retention <- ggplot(retention,
                          aes(x = MeanCoverage, y = prop, colour = method)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    x = "Mean sample coverage (x)",
    y = "Proportion of ANGSD GL hets retained",
    colour = "Variant caller",
    title = "Does sample coverage explain retained ANGSD GL heterozygous sites?"
  ) +
  theme(legend.position = "bottom")

print(p_cov_retention)

ggsave(file.path(outdir, "06_retained_het_vs_sample_coverage.png"),
       p_cov_retention, width = 10, height = 6, dpi = 300)


# ------------------------------------------------------------
# 15. Supervisor summary figure
# ------------------------------------------------------------

combined <- combined <- (p1 / p_depth / p_balance) +
  plot_annotation(
    title = "Diagnostic comparison of ANGSD GL heterozygous sites across VCF callers"
  )

ggsave(file.path(outdir, "supervisor_summary_figure_all_callers.png"),
       combined, width = 11, height = 13, dpi = 300)

# ------------------------------------------------------------
# 16. Console summaries
# ------------------------------------------------------------

message("Overall outcome by caller:")

summary_by_method %>%
  arrange(method, desc(prop)) %>%
  print(n = Inf)

message("Diagnostics by caller and outcome:")

summary_mechanism %>%
  arrange(method, outcome_clean) %>%
  print(n = Inf)