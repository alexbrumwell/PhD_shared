# SCRIPT TEST HET SCAF1 ANGSD VS SNPAD

rm(list = ls())
# To view plots 
library(httpgd)
hgd()
# At the end type: dev.off()
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(writexl)
library(readxl)


setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Het/tests_scaf1_angsdvssnpad")

df <- read_xlsx("summary_sites_scaf1.xlsx")

# Add metadata 
metadata <- read_excel(
  "~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/OBP_sumstats.xlsx",
  sheet = 1
) %>%
  rename(Individual = Sample_ID) 

df <- df %>% left_join(metadata, by = "Individual")

# Clean numeric columns
num_cols <- c(
  "Raw_Scaf1_Count", "Luke_Scaf1_Count", "ANGSD_Raw_Scaf1_Count",
  "ANGSD_Total_Sites", "ANGSD_Het_Sites",
  "SNPAD_Raw_Scaf1_Count", "SNPAD_Total_Sites", "SNPAD_Het_Sites",
  "ANGSD_Total_Sites_maf0.02", "ANGSD_Het_Sites_maf0.02", "Het_maf0.02",
  "ANGSD_Total_Sites_maf0.05", "ANGSD_Het_Sites_maf0.05", "Het_maf0.05",
  "MeanCoverage"
)

df <- df %>%
  mutate(across(any_of(num_cols), as.numeric))

# Derived columns
df2 <- df %>%
  mutate(
    ANGSD_het_rate = ANGSD_Het_Sites / ANGSD_Total_Sites,
    SNPAD_het_rate = SNPAD_Het_Sites / SNPAD_Total_Sites,

    Het_diff_n = ANGSD_Het_Sites - SNPAD_Het_Sites,
    Het_diff_pct = if_else(
      ANGSD_Het_Sites > 0,
      100 * (SNPAD_Het_Sites - ANGSD_Het_Sites) / ANGSD_Het_Sites,
      NA_real_
    ),
    Site_diff_n = ANGSD_Total_Sites - SNPAD_Total_Sites,
    Site_retention_snpAD_vs_ANGSD = if_else(
      ANGSD_Total_Sites > 0,
      SNPAD_Total_Sites / ANGSD_Total_Sites,
      NA_real_
    ),
    Het_retention_snpAD_vs_ANGSD = if_else(
      ANGSD_Het_Sites > 0,
      SNPAD_Het_Sites / ANGSD_Het_Sites,
      NA_real_
    ),

    Het_change_maf0.02_n = ANGSD_Het_Sites_maf0.02 - ANGSD_Het_Sites,
    Het_change_maf0.05_n = ANGSD_Het_Sites_maf0.05 - ANGSD_Het_Sites,

    Het_change_maf0.02_pct = if_else(
      ANGSD_Het_Sites > 0,
      100 * (ANGSD_Het_Sites_maf0.02 - ANGSD_Het_Sites) / ANGSD_Het_Sites,
      NA_real_
    ),
    Het_change_maf0.05_pct = if_else(
      ANGSD_Het_Sites > 0,
      100 * (ANGSD_Het_Sites_maf0.05 - ANGSD_Het_Sites) / ANGSD_Het_Sites,
      NA_real_
    )
  )

# Order samples by year if a year column exists
if ("Year" %in% names(df2)) {
  sample_order <- df2 %>%
    distinct(Individual, Group, Year) %>%
    arrange(Group, Year, Individual) %>%
    pull(Individual)
} 
df2 <- df2 %>%
  mutate(Individual = factor(Individual, levels = sample_order))

# Save full summary table
write_xlsx(df2, "summary_sites_scaf1_with_metrics.xlsx")


# Long format for plotting
plot_sites <- df2 %>%
  select(Individual, Group, MeanCoverage, ANGSD_Total_Sites, SNPAD_Total_Sites) %>%
  pivot_longer(
    cols = c(ANGSD_Total_Sites, SNPAD_Total_Sites),
    names_to = "Method",
    values_to = "Total_Sites"
  )

plot_hets <- df2 %>%
  select(Individual, Group, MeanCoverage, ANGSD_Het_Sites, SNPAD_Het_Sites) %>%
  pivot_longer(
    cols = c(ANGSD_Het_Sites, SNPAD_Het_Sites),
    names_to = "Method",
    values_to = "Het_Sites"
  )

# Long format for ANGSD across AF filters
plot_maf_hets <- df2 %>%
  select(
    Individual, Group, MeanCoverage,
    ANGSD_Het_Sites,
    ANGSD_Het_Sites_maf0.02,
    ANGSD_Het_Sites_maf0.05
  ) %>%
  pivot_longer(
    cols = c(
      ANGSD_Het_Sites,
      ANGSD_Het_Sites_maf0.02,
      ANGSD_Het_Sites_maf0.05
    ),
    names_to = "Filter",
    values_to = "Het_Sites"
  ) %>%
  mutate(
    Filter = recode(
      Filter,
      "ANGSD_Het_Sites" = "No filter",
      "ANGSD_Het_Sites_maf0.02" = "MAF 0.02",
      "ANGSD_Het_Sites_maf0.05" = "MAF 0.05"
    ),
    Filter = factor(Filter, levels = c("No filter", "MAF 0.02", "MAF 0.05"))
  )

plot_maf_rates <- df2 %>%
  select(
    Individual, Group, MeanCoverage,
    ANGSD_het_rate,
    Het_maf0.02,
    Het_maf0.05
  ) %>%
  pivot_longer(
    cols = c(ANGSD_het_rate, Het_maf0.02, Het_maf0.05),
    names_to = "Filter",
    values_to = "Het_Rate"
  ) %>%
  mutate(
    Filter = recode(
      Filter,
      "ANGSD_het_rate" = "No filter",
      "Het_maf0.02" = "MAF 0.02",
      "Het_maf0.05" = "MAF 0.05"
    ),
    Filter = factor(Filter, levels = c("No filter", "MAF 0.02", "MAF 0.05"))
  )


# Coverage label data for total-sites plot
df_labels_sites <- plot_sites %>%
  group_by(Individual, Group) %>%
  summarise(
    y_pos = max(Total_Sites, na.rm = TRUE) * 1.05,
    MeanCoverage = first(MeanCoverage),
    .groups = "drop"
  ) %>%
  mutate(
    sample_label = if_else(
      MeanCoverage < 4,
      paste0("<span style='color:red;'>", Individual, "</span>"),
      as.character(Individual)
    )
  )

# Coverage label data for het-sites plot
df_labels_hets <- plot_hets %>%
  group_by(Individual, Group) %>%
  summarise(
    y_pos = max(Het_Sites, na.rm = TRUE) * 1.05,
    MeanCoverage = first(MeanCoverage),
    .groups = "drop"
  ) %>%
  mutate(
    sample_label = if_else(
      MeanCoverage < 4,
      paste0("<span style='color:red;'>", Individual, "</span>"),
      as.character(Individual)
    )
  )

# Coverage labels for MAF het-sites plot
df_labels_maf_hets <- plot_maf_hets %>%
  group_by(Individual, Group) %>%
  summarise(
    y_pos = max(Het_Sites, na.rm = TRUE) * 1.05,
    MeanCoverage = first(MeanCoverage),
    .groups = "drop"
  )



# Build coloured x-axis labels
label_df <- df2 %>%
  distinct(Individual, MeanCoverage) %>%
  mutate(
    axis_label = if_else(
      MeanCoverage < 4,
      paste0("<span style='color:red;'>", Individual, "</span>"),
      as.character(Individual)
    )
  )

label_vec <- setNames(label_df$axis_label, label_df$Individual)

# Plot A: total sites
p1 <- ggplot(plot_sites, aes(x = Individual, y = Total_Sites, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(
    data = df_labels_sites,
    aes(x = Individual, y = y_pos, label = sprintf("%.1f×", MeanCoverage)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Scaffold 1: Total sites by method using same list of sites (Luke's list)",
    x = "Sample",
    y = "Total analysed sites",
    subtitle = "Samples are ordered chronologically. Coverage shown above bars. Samples with coverage <4× are coloured red."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(p1)
ggsave("plot_total_sites_scaf1.png", p1, width = 12, height = 8, dpi = 300)


# Plot B: heterozygous sites
p2 <- ggplot(plot_hets, aes(x = Individual, y = Het_Sites, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(
    data = df_labels_hets,
    aes(x = Individual, y = y_pos, label = sprintf("%.1f×", MeanCoverage)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Scaffold 1: Heterozygous sites by method using same list of sites (Luke's list)",
    x = "Sample",
    y = "Heterozygous sites",
    subtitle = "amples are ordered chronologically. Coverage shown above bars. Samples with coverage <4× are coloured red."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(p2)
ggsave("plot_het_sites_scaf1.png", p2, width = 12, height = 8, dpi = 300)


# Plot C: retention of sites vs retention of hets
p3 <- ggplot(df2, aes(x = Site_retention_snpAD_vs_ANGSD, y = Het_retention_snpAD_vs_ANGSD)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    title = "Does snpAD lose hets more than expected from site loss?",
    x = "snpAD / ANGSD total-site retention",
    y = "snpAD / ANGSD het-site retention"
  )

print(p3)
ggsave("plot_retention_scaf1.png", p3, width = 6, height = 5, dpi = 300)


# PLOT BOXPLOT OF SITES
plot_box_combined <- df2 %>%
  select(
    Individual, Group, MeanCoverage,
    ANGSD_Total_Sites, SNPAD_Total_Sites,
    ANGSD_Het_Sites, SNPAD_Het_Sites
  ) %>%
  pivot_longer(
    cols = c(
      ANGSD_Total_Sites, SNPAD_Total_Sites,
      ANGSD_Het_Sites, SNPAD_Het_Sites
    ),
    names_to = c("Method", "Metric"),
    names_pattern = "(ANGSD|SNPAD)_(.*)"
  ) %>%
  mutate(
    Method = recode(Method, ANGSD = "ANGSD", SNPAD = "snpAD"),
    Metric = recode(
      Metric,
      Total_Sites = "Total analysed sites",
      Het_Sites = "Heterozygous sites"
    )
  )


  )

print(p_box_grid)
ggsave("plot_box_grid_sites_hets_scaf1.png", p_box_grid, width = 10, height = 7, dpi = 300)

# PLOT D: does the difference in het sites between methods depend on coverage?
# PLOT: % difference in het sites vs coverage
p_diff_pct <- ggplot(
  df2,
  aes(
    x = MeanCoverage,
    y = Het_diff_pct,
    colour = Group
  )
) +
  geom_point(
    aes(shape = MeanCoverage < 4),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("≥4×", "<4×"),
    name = "Coverage"
  ) +
  labs(
  title = "Percent difference in heterozygous sites between methods vs coverage",
  subtitle = "Negative values indicate fewer het sites in snpAD compared to ANGSD",
  x = "Mean coverage (×)",
  y = "% diff in het sites of snpAD compared to ANGSD"
) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

print(p_diff_pct)
ggsave("plot_coverage_vs_het_diff_pct_scaf1.png", p_diff_pct, width = 7, height = 5, dpi = 300)


# PLOT: HOW DOES COVERAGE AFFECT HET BY METHOD?
plot_het_rate <- df2 %>%
  select(
    Individual, Group, MeanCoverage,
    ANGSD_het_rate, SNPAD_het_rate
  ) %>%
  pivot_longer(
    cols = c(ANGSD_het_rate, SNPAD_het_rate),
    names_to = "Method",
    values_to = "het_rate"
  ) %>%
  mutate(
    Method = recode(
      Method,
      ANGSD_het_rate = "ANGSD",
      SNPAD_het_rate = "snpAD"
    ),
    CoverageClass = if_else(MeanCoverage < 4, "<4×", "≥4×")
  )

p4b <- ggplot(
  plot_het_rate,
  aes(
    x = MeanCoverage,
    y = het_rate,
    colour = Method
  )
) +
  geom_point(
    aes(shape = CoverageClass),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(~ Group, scales = "free_x") +
  scale_shape_manual(
    values = c("≥4×" = 16, "<4×" = 17),
    name = "Coverage"
  ) +
  labs(
    title = "Effect of coverage on heterozygosity by method",
    subtitle = "Scaffold 1, using Luke's site list. Heterozygosity is calculated as heterozygous sites / total analysed sites.",
    x = "Mean coverage (×)",
    y = "Heterozygosity rate"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

print(p4b)
ggsave("plot_coverage_vs_het_rate_by_group_scaf1.png", p4b, width = 8, height = 5, dpi = 300)


# PLOT: MAF FILTERS EFFECT
p_maf_hets <- ggplot(plot_maf_hets, aes(x = Individual, y = Het_Sites, fill = Filter)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(
    data = df_labels_maf_hets,
    aes(x = Individual, y = y_pos, label = sprintf("%.1f×", MeanCoverage)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Scaffold 1: Heterozygous sites across allele-frequency filters",
    subtitle = "ANGSD results only. Coverage shown above bars. Samples with coverage <4× are coloured red.",
    x = "Sample",
    y = "Heterozygous sites",
    fill = "AF filter"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(p_maf_hets)
ggsave("plot_het_sites_by_af_filter_scaf1.png", p_maf_hets, width = 13, height = 8, dpi = 300)


p_box_maf_hets <- ggplot(
  plot_maf_hets,
  aes(x = Filter, y = Het_Sites, colour = Group)
) +
  geom_boxplot(
    aes(group = interaction(Filter, Group)),
    width = 0.7,
    alpha = 0.8,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    colour = "black"
  ) +
  geom_point(
    aes(
      shape = MeanCoverage < 4,
      group = interaction(Filter, Group)
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2,
    alpha = 0.7
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("≥4×", "<4×"),
    name = "Coverage"
  ) +
  facet_wrap(~ Group, ) +
  labs(
    title = "Change in heterozygous sites across allele-frequency filters",
    subtitle = "ANGSD results only, scaffold 1",
    x = "",
    y = "Heterozygous sites"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_box_maf_hets)
ggsave("plot_box_het_sites_by_af_filter_scaf1.png", p_box_maf_hets, width = 9, height = 5, dpi = 300)

# Coverage vs het sites
p_cov_vs_hets_facet <- ggplot(
  plot_maf_hets,
  aes(
    x = MeanCoverage,
    y = Het_Sites
  )
) +
  geom_point(
    aes(
      colour = Group,
      shape = MeanCoverage < 4
    ),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_smooth(
    aes(colour = Group),
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_grid(Filter ~ Group,scales = "free_x") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("≥4×", "<4×"),
    name = "Coverage"
  ) +
  labs(
    title = "Effect of coverage on heterozygous site counts across MAF filters",
    subtitle = "Rows = MAF filter, Columns = Group (ANGSD results, scaffold 1)",
    x = "Mean coverage (×)",
    y = "Heterozygous sites",
    colour = "Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

print(p_cov_vs_hets_facet)

ggsave(
  "plot_coverage_vs_het_sites_maf_faceted_rows_filter_cols_group.png",
  p_cov_vs_hets_facet,
  width = 12,
  height = 10,
  dpi = 300
)


# PLOT: MAF EFFECTS INCLUDING SNPAD
# Long format including snpAD and ANGSD AF-filtered results
plot_hets_all <- df2 %>%
  select(
    Individual, Group, MeanCoverage,
    SNPAD_Het_Sites,
    ANGSD_Het_Sites,
    ANGSD_Het_Sites_maf0.02,
    ANGSD_Het_Sites_maf0.05
  ) %>%
  pivot_longer(
    cols = c(
      SNPAD_Het_Sites,
      ANGSD_Het_Sites,
      ANGSD_Het_Sites_maf0.02,
      ANGSD_Het_Sites_maf0.05
    ),
    names_to = "Filter",
    values_to = "Het_Sites"
  ) %>%
  mutate(
    Filter = recode(
      Filter,
      "SNPAD_Het_Sites" = "snpAD",
      "ANGSD_Het_Sites" = "ANGSD",
      "ANGSD_Het_Sites_maf0.02" = "ANGSD + MAF 0.02",
      "ANGSD_Het_Sites_maf0.05" = "ANGSD + MAF 0.05"
    ),
    Filter = factor(
      Filter,
      levels = c("ANGSD", "ANGSD + MAF 0.02", "ANGSD + MAF 0.05", "snpAD")
    )
  )


df_labels_hets_all <- plot_hets_all %>%
  group_by(Individual, Group) %>%
  summarise(
    y_pos = max(Het_Sites, na.rm = TRUE) * 1.05,
    MeanCoverage = first(MeanCoverage),
    .groups = "drop"
  )


# PLOT: snpAD plus ANGSD AF filters
p_hets_all <- ggplot(plot_hets_all, aes(x = Individual, y = Het_Sites, fill = Filter)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.65) +
  geom_text(
    data = df_labels_hets_all,
    aes(x = Individual, y = y_pos, label = sprintf("%.1f×", MeanCoverage)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Scaffold 1: Heterozygous sites across methods and allele-frequency filters",
    subtitle = "Coverage shown above bars. Samples with coverage <4× are coloured red.",
    x = "Sample",
    y = "Heterozygous sites",
    fill = "Method / filter"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(p_hets_all)
ggsave("plot_het_sites_methods_and_af_filters_scaf1.png", p_hets_all, width = 14, height = 8, dpi = 300)


# BOXPLOT: snpAD plus ANGSD AF filters
p_box_hets_all <- ggplot(
  plot_hets_all,
  aes(x = Filter, y = Het_Sites, colour = Group)
) +
  geom_boxplot(
    aes(group = interaction(Filter, Group)),
    width = 0.7,
    alpha = 0.8,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    colour = "black"
  ) +
  geom_point(
    aes(
      shape = MeanCoverage < 4,
      group = interaction(Filter, Group)
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2,
    alpha = 0.7
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("≥4×", "<4×"),
    name = "Coverage"
  ) +
  facet_wrap(~ Group, ) +
  labs(
    title = "Comparison of heterozygous sites across methods and allele-frequency filters",
    subtitle = "Scaffold 1",
    x = "",
    y = "Heterozygous sites"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95")
  )

print(p_box_hets_all)
ggsave(
  "plot_box_het_sites_methods_and_af_filters_scaf1.png",
  p_box_hets_all,
  width = 10,
  height = 8,
  dpi = 300
)


# PLOT FOR CLAUDIA: ANGSD and SNPAD with two denominator types

library(ggrepel)

df2 <- df2 %>%
  mutate(
    ANGSD_rate_individual_sites = ANGSD_Het_Sites / ANGSD_Total_Sites,
    ANGSD_rate_luke_sites       = ANGSD_Het_Sites / Luke_Scaf1_Count,
    SNPAD_rate_individual_sites = SNPAD_Het_Sites / SNPAD_Total_Sites,
    SNPAD_rate_luke_sites       = SNPAD_Het_Sites / Luke_Scaf1_Count,
    ratio_ANGSD_snpAD_heterozygosity = ANGSD_het_rate / SNPAD_het_rate,
    dif_ANGSD_snpAD_het = ANGSD_het_rate - SNPAD_het_rate,
    samesite_ratio_ANGSD_snpAD_heterozygosity = ANGSD_rate_luke_sites / SNPAD_rate_luke_sites,
    samesite_dif_ANGSD_snpAD_het = ANGSD_rate_luke_sites - SNPAD_rate_luke_sites
    
  )

# Order individuals by coverage within group
sample_order_cov <- df2 %>%
  arrange(Group, MeanCoverage, Individual) %>%
  pull(Individual)

df2 <- df2 %>%
  mutate(
    Individual = factor(Individual, levels = sample_order_cov)
  )

plot_het_compare <- df2 %>%
  select(
    Individual,
    Group,
    MeanCoverage,
    ANGSD_rate_individual_sites,
    ANGSD_rate_luke_sites,
    SNPAD_rate_individual_sites,
    SNPAD_rate_luke_sites
  ) %>%
  pivot_longer(
    cols = c(
      ANGSD_rate_individual_sites,
      ANGSD_rate_luke_sites,
      SNPAD_rate_individual_sites,
      SNPAD_rate_luke_sites
    ),
    names_to = "Estimate_type",
    values_to = "Het_Rate"
  ) %>%
  separate(
    Estimate_type,
    into = c("Method", "Denominator"),
    sep = "_rate_"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("ANGSD", "SNPAD")),
    Denominator = recode(
      Denominator,
      "individual_sites" = "Per individual callable sites",
      "luke_sites" = "Total callable sites (same across samples)"
    ),
    Denominator = factor(
      Denominator,
      levels = c(
        "Per individual callable sites",
        "Total callable sites (same across samples)"
      )
    )
  )

df_labels_compare <- plot_het_compare %>%
  group_by(Individual, Group, Method) %>%
  summarise(
    y_pos = max(Het_Rate, na.rm = TRUE) * 1.05,
    MeanCoverage = first(MeanCoverage),
    .groups = "drop"
  )

p_compare <- ggplot(
  plot_het_compare,
  aes(
    x = Individual,
    y = Het_Rate,
    fill = Denominator
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    data = df_labels_compare,
    aes(
      x = Individual,
      y = y_pos,
      label = sprintf("%.1f×", MeanCoverage)
    ),
    inherit.aes = FALSE,
    size = 2.6,
    vjust = 0
  ) +
  facet_grid(Method ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Comparison of ANGSD and SNPAD heterozygosity estimates",
    subtitle = "Individuals ordered by increasing coverage. Each method is shown with two denominator choices.",
    x = "Individual",
    y = "Heterozygosity estimate",
    fill = "Denominator"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(p_compare)

ggsave(
  "Claudia_plot_comparison_ANGSD_SNPAD_het_estimates_scaf1.png",
  p_compare,
  width = 14,
  height = 9,
  dpi = 300
)

p_cov_vs_het <- ggplot(
  plot_het_compare,
  aes(
    x = MeanCoverage,
    y = Het_Rate,
    colour = Group
  )
) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_grid(Method ~ Denominator, scales = "fixed") +
  labs(
    title = "Does coverage influence ANGSD and SNPAD heterozygosity estimates?",
    subtitle = "Each method is shown using either per-individual callable sites or a shared Luke scaffold 1 denominator",
    x = "Mean coverage",
    y = "Heterozygosity estimate",
    colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom"
  )

p_cov_vs_het_labelled <- p_cov_vs_het +
  ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 3,
    show.legend = FALSE,
    max.overlaps = 1
  )

print(p_cov_vs_het_labelled)

ggsave(
  "Claudia_plot_cov_vs_ANGSD_SNPAD_het_estimates_scaf1.png",
  p_cov_vs_het_labelled,
  width = 14,
  height = 9,
  dpi = 300
)
    

ggplot(df2, aes(x = Individual, y = ratio_ANGSD_snpAD_heterozygosity)) +
  geom_col(position = "dodge")



p_claudia_1 <- ggplot(df2, aes(x = Individual, y = dif_ANGSD_snpAD_het)) +
  geom_col(position = "dodge")
p_claudia_1 <- ggsave("tmp_Claudia_fig1.png", p_claudia_1)


ggplot(data = df2 %>% filter(MeanCoverage > 4), aes(x = Individual, y = samesite_ratio_ANGSD_snpAD_heterozygosity, fill = MeanCoverage) )+
  geom_col(position = "dodge") 


ggplot(df2, aes(x = Individual, y = samesite_dif_ANGSD_snpAD_het)) +
  geom_col(position = "dodge") 


