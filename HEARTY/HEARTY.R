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

