# SCRIPT TEST HET SCAF1 ANGSD VS SNPAD
# KEY QUESTION: Does downsampling to 5× affect heterozygosity estimates differently for ANGSD vs SNPAD?

rm(list = ls())

library(httpgd)
hgd()

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(writexl)
library(readxl)
library(stringr)
library(patchwork)


setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Het/tests_scaf1_angsdvssnpad/down5x")

# Read in data 
df_down <- read_xlsx("summary_sites_scaf1.xlsx") %>%
  rename_with(~ paste0(., "_down"), -Individual)

df_orig <- read_xlsx("../summary_sites_scaf1.xlsx") %>%
  rename_with(~ paste0(., "_orig"), -Individual)

# Join original and downsampled data
df_joined <- df_down %>%
  inner_join(df_orig, by = "Individual")

# Calculate differences
df_compare <- df_joined %>%
  mutate(
    ANGSD_Het_sites_diff = ANGSD_Het_Sites_down - ANGSD_Het_Sites_orig,
    SNPAD_Het_sites_diff = SNPAD_Het_Sites_down - SNPAD_Het_Sites_orig,

    ANGSD_Het_sites_ratio = ANGSD_Het_Sites_down / ANGSD_Het_Sites_orig,
    SNPAD_Het_sites_ratio = SNPAD_Het_Sites_down / SNPAD_Het_Sites_orig,

    ANGSD_site_retention = ANGSD_Total_Sites_down / ANGSD_Total_Sites_orig,
    SNPAD_site_retention = SNPAD_Total_Sites_down / SNPAD_Total_Sites_orig
  )

# Add metadata 
metadata <- read_excel(
  "~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/OBP_sumstats.xlsx",
  sheet = 1
) %>%
  rename(Individual = Sample_ID)

df <- df_compare %>%
  left_join(metadata, by = "Individual") 

df$MeanCoverage_down5x <- as.numeric(df$MeanCoverage_down5x)

# Remove any individuals that were NOT downsampled
df_filt <- df %>%
  filter(!is.na(MeanCoverage_down5x)) %>%
  mutate(
    Year_num = suppressWarnings(as.numeric(Year))
  )





# QUESTION 1: IS DOWNSAMPLED HET MORE RELIABLE?
# Build long dataframe with het estimate
plot_het_bar <- bind_rows(
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year_num,
      Treatment = "Original",
      Method = "ANGSD",
      Het_Sites = ANGSD_Het_Sites_orig,
      Total_Sites = ANGSD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year_num,
      Treatment = "Downsampled",
      Method = "ANGSD",
      Het_Sites = ANGSD_Het_Sites_down,
      Total_Sites = ANGSD_Total_Sites_down
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year_num,
      Treatment = "Original",
      Method = "SNPAD",
      Het_Sites = SNPAD_Het_Sites_orig,
      Total_Sites = SNPAD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year_num,
      Treatment = "Downsampled",
      Method = "SNPAD",
      Het_Sites = SNPAD_Het_Sites_down,
      Total_Sites = SNPAD_Total_Sites_down
    )
) %>%
  mutate(
    Het_Estimate = Het_Sites / Total_Sites,
    Treatment = factor(Treatment, levels = c("Original", "Downsampled")),
    Method = factor(Method, levels = c("ANGSD", "SNPAD")),
    Method_Treatment = factor(
      paste(Method, Treatment),
      levels = c(
        "ANGSD Original",
        "ANGSD Downsampled",
        "SNPAD Original",
        "SNPAD Downsampled"
      )
    )
  )

# Order individuals by group, year, then ID
sample_order <- df_filt %>%
  arrange(Group, Year_num, Individual) %>%
  pull(Individual)

plot_het_bar$Individual <- factor(plot_het_bar$Individual, levels = sample_order)

df_labels_cov <- plot_het_bar %>%
  group_by(Individual, Group) %>%
  summarise(
    y_pos = max(Het_Estimate, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  left_join(
    df_filt %>%
      select(Individual, MeanCoverage),
    by = "Individual"
  ) %>%
  mutate(
cov_label = paste0(formatC(MeanCoverage, format = "f", digits = 1), "×")  )

# Create red coloured indivudal labels for individuals that have higher het in the downsampled treatementment compared to original
label_df <- df_filt %>%
  mutate(
    ANGSD_outlier =
      (ANGSD_Het_Sites_down / ANGSD_Total_Sites_down) >
      (ANGSD_Het_Sites_orig / ANGSD_Total_Sites_orig) &
      ANGSD_Het_Sites_down < ANGSD_Het_Sites_orig,
    label = if_else(
      ANGSD_outlier,
      paste0("<span style='color:red;'>", Individual, "</span>"),
      Individual
    )
  )

label_vec <- setNames(label_df$label, label_df$Individual)

# ----------------------------
# Plot
# ----------------------------
p_het_bar_cov <- ggplot(
  plot_het_bar,
  aes(x = Individual, y = Het_Estimate, fill = Method_Treatment)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    data = df_labels_cov,
    aes(x = Individual, y = y_pos, label = cov_label),
    inherit.aes = FALSE,
    size = 2.7,
    vjust = 0,
    lineheight = 0.9
  ) +
  facet_grid(Method ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Heterozygosity estimates by method and downsampling to 5x",
    subtitle = "Samples that have not been downsampled are not shown.Original coverage shown above each sample",
    x = "Sample",
    y = "Heterozygosity estimate",
    fill = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

print(p_het_bar_cov)

p_het_count_bar_cov <- ggplot(
  plot_het_bar,
  aes(x = Individual, y = Het_Sites, fill = Method_Treatment)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    data = df_labels_cov,
    aes(x = Individual, y = y_pos, label = cov_label),
    inherit.aes = FALSE,
    size = 2.7,
    vjust = 0,
    lineheight = 0.9
  ) +
  facet_grid(Method ~ Group, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = label_vec) +
  labs(
    title = "Heterozygous site counts by method and downsampling to 5x",
    subtitle = "Samples that have not been downsampled are not shown.Original coverage shown above each sample",
    x = "Sample",
    y = "Heterozygous site count",
    fill = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )
print(p_het_count_bar_cov)

p_combined <- (p_het_bar_cov + p_het_count_bar_cov ) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("plot_bar_heterozygosity_down5x_combined.png", p_combined, width = 20, height = 10, dpi = 300)


# PLOT HET SITES VS HET ESTIMATE
plot_scatter_df <- bind_rows(
  df_filt %>%
    transmute(
      Individual,
      Group,
      Method = "ANGSD",
      Treatment = "Original",
      Het_Sites = ANGSD_Het_Sites_orig,
      Total_Sites = ANGSD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Method = "ANGSD",
      Treatment = "Downsampled",
      Het_Sites = ANGSD_Het_Sites_down,
      Total_Sites = ANGSD_Total_Sites_down
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Method = "SNPAD",
      Treatment = "Original",
      Het_Sites = SNPAD_Het_Sites_orig,
      Total_Sites = SNPAD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Method = "SNPAD",
      Treatment = "Downsampled",
      Het_Sites = SNPAD_Het_Sites_down,
      Total_Sites = SNPAD_Total_Sites_down
    )
) %>%
  mutate(
    Het_Estimate = Het_Sites / Total_Sites,
    Treatment = factor(Treatment, levels = c("Original", "Downsampled")),
    Method = factor(Method, levels = c("ANGSD", "SNPAD"))
  )

mean_totals <- plot_scatter_df %>%
  group_by(Method, Treatment) %>%
  summarise(mean_total = mean(Total_Sites), .groups = "drop")

p_scatter <- ggplot(
  plot_scatter_df,
  aes(x = Het_Sites, y = Het_Estimate)
) +
  geom_point(aes(colour = Group), size = 2.5, alpha = 0.8) +
  geom_abline(
    data = mean_totals,
    aes(intercept = 0, slope = 1 / mean_total),
    linetype = 2,
    colour = "grey40"
  ) +
  ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 4,
    show.legend = FALSE
  ) +
  facet_grid(Method ~ Treatment) +
  labs(
    x = "Heterozygous sites",
    y = "Heterozygosity estimate",
    title = "Het sites vs heterozygosity estimate",
    subtitle = "Dashed line = expected relationship given mean total sites per group and treatment"
  ) +
  theme_bw()

print(p_scatter)
ggsave("plot_scatter_het_sites_vs_estimate.png", p_scatter, width = 12, height = 8, dpi = 300)
 
# PLOT HET SITES VS TOTAL SITES
df_ratio <- df_filt %>%
  transmute(
    Individual,
    Group,

    ANGSD_het_ret = ANGSD_Het_Sites_down / ANGSD_Het_Sites_orig,
    ANGSD_total_ret = ANGSD_Total_Sites_down / ANGSD_Total_Sites_orig,

    SNPAD_het_ret = SNPAD_Het_Sites_down / SNPAD_Het_Sites_orig,
    SNPAD_total_ret = SNPAD_Total_Sites_down / SNPAD_Total_Sites_orig
  ) %>%
  pivot_longer(
    cols = -c(Individual, Group),
    names_to = "var",
    values_to = "value"
  ) %>%
  separate(var, into = c("Method", "Metric"), sep = "_") %>%
  pivot_wider(names_from = Metric, values_from = value)

p_het_total <- ggplot(df_ratio, aes(x = total, y = het, colour = Group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(colour = Group), size = 2.5) +
    ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 3,
    show.legend = FALSE
  ) +
  labs(
    x = "Total site retention (downsampled / original)",
    y = "Het site retention (downsampled / original)",
    title = "Relative loss of total vs heterozygous sites after downsampling to 5x",
    subtitle = "Reference line is 1:1. Points above dashed line = het sites are retained at a higher rate compared to total sites after downsampling"
  ) +
  facet_grid(. ~ Method) +
  theme_bw()
print(p_het_total)
ggsave("plot_het_vs_total_retention.png", p_het_total, width = 10, height = 6, dpi = 300)

  
# ----------------------------
# Build long dataframe for plotting
# ----------------------------
plot_df <- bind_rows(
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year,
      MeanCoverage,
      Treatment = "Original",
      Method = "ANGSD",
      Het_Sites = ANGSD_Het_Sites_orig,
      Total_Sites = ANGSD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year,
      MeanCoverage,
      Treatment = "Original",
      Method = "SNPAD",
      Het_Sites = SNPAD_Het_Sites_orig,
      Total_Sites = SNPAD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year,
      MeanCoverage,
      Treatment = "Downsampled",
      Method = "ANGSD",
      Het_Sites = ANGSD_Het_Sites_down,
      Total_Sites = ANGSD_Total_Sites_down
    ),
  df_filt %>%
    transmute(
      Individual,
      Group,
      Year,
      MeanCoverage,
      Treatment = "Downsampled",
      Method = "SNPAD",
      Het_Sites = SNPAD_Het_Sites_down,
      Total_Sites = SNPAD_Total_Sites_down
    )
) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Original", "Downsampled")),
    Method = factor(Method, levels = c("ANGSD", "SNPAD")),
    CoverageClass = if_else(MeanCoverage < 4, "<4×", "≥4×"),
    Het_Rate = Het_Sites / Total_Sites
  )


# QUESTION 2: DOES DOWNSAMPLING CHANGE HET SITES?
plot_change_df <- bind_rows(
  df_filt %>%
    transmute(
      Individual, Group, Method = "ANGSD", Treatment = "Original",
      Het_Sites = ANGSD_Het_Sites_orig,
      Het_Estimate = ANGSD_Het_Sites_orig / ANGSD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual, Group, Method = "ANGSD", Treatment = "Downsampled",
      Het_Sites = ANGSD_Het_Sites_down,
      Het_Estimate = ANGSD_Het_Sites_down / ANGSD_Total_Sites_down
    ),
  df_filt %>%
    transmute(
      Individual, Group, Method = "SNPAD", Treatment = "Original",
      Het_Sites = SNPAD_Het_Sites_orig,
      Het_Estimate = SNPAD_Het_Sites_orig / SNPAD_Total_Sites_orig
    ),
  df_filt %>%
    transmute(
      Individual, Group, Method = "SNPAD", Treatment = "Downsampled",
      Het_Sites = SNPAD_Het_Sites_down,
      Het_Estimate = SNPAD_Het_Sites_down / SNPAD_Total_Sites_down
    )
) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Original", "Downsampled")),
    Method = factor(Method, levels = c("ANGSD", "SNPAD"))
  )
increase_df <- plot_change_df %>%
  select(Individual, Method, Treatment, Het_Estimate) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = Het_Estimate
  ) %>%
  mutate(
    Increased = Downsampled > Original
  ) %>%
  select(Individual, Method, Increased)
plot_change_df2 <- plot_change_df %>%
  left_join(increase_df, by = c("Individual", "Method"))
# Build long dataframe with both metrics
plot_change_long <- plot_change_df2 %>%
  pivot_longer(
    cols = c(Het_Sites, Het_Estimate),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      Het_Sites = "Heterozygous sites",
      Het_Estimate = "Heterozygosity estimate"
    ),
    Metric = factor(
      Metric,
      levels = c("Heterozygous sites", "Heterozygosity estimate")
    )
  )

p_combined <- ggplot(
  plot_change_long,
  aes(x = Treatment, y = Value, group = interaction(Individual, Method))
) +
  geom_line(
    aes(colour = ifelse(Metric == "Heterozygosity estimate", Increased, NA)),
    alpha = 0.6,
    linewidth = 0.4,
    show.legend = TRUE
  ) +
  geom_point(
    aes(colour = Group),
    size = 2,
  ) +
  facet_grid(Metric ~ Method, scales = "free_y" ) +
  # Two separate colour scales merged manually
  scale_colour_manual(
    values = c(
      "TRUE" = "red",
      "FALSE" = "grey80",
      "Historical" = "#F8766D",
      "Modern" = "#00BFC4"
    ),
    breaks = c("Historical", "Modern", "TRUE"),
    labels = c(
      "Historical",
      "Modern",
      "Increase in het estimate"
    ),
    name = ""
  ) +
  labs(
    x = "",
    y = "",
    title = "Effect of downsampling on heterozygous sites and heterozygosity estimates",
    subtitle = "Red lines indicate individuals where heterozygosity estimate increased after downsampling"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

print(p_combined)

ggsave("plot_down5x_het_sites_and_estimate_faceted.png", p_combined, width = 10, height = 8)

# QUESTION 3: DOES METHOD HAVE EFFECT ON DOWNSAMPLING? DOES METHOD AFFECT TOAL SITES AS WELL AS HET SITES?
# RESULT: Downsampling primarily affects heterozygosity estimation
# Create retention rate of downsampled het sites compared to original, and compare to retention of total sites to see if het sites are disproportionately affected
ret_df <- df_filt %>%
  transmute(
    Individual,
    Group,
    MeanCoverage_orig = MeanCoverage,
    MeanCoverage_down = MeanCoverage_down5x,

    ANGSD = ANGSD_Het_Sites_down / ANGSD_Het_Sites_orig,
    SNPAD = SNPAD_Het_Sites_down / SNPAD_Het_Sites_orig
  ) %>%
  pivot_longer(
    cols = c(ANGSD, SNPAD),
    names_to = "Method",
    values_to = "Het_retention"
  )
ret_df2 <- df_filt %>%
  transmute(
    Individual,
    Group,
    MeanCoverage_down = MeanCoverage_down5x,

    ANGSD_Het = ANGSD_Het_Sites_down / ANGSD_Het_Sites_orig,
    SNPAD_Het = SNPAD_Het_Sites_down / SNPAD_Het_Sites_orig,

    ANGSD_Total = ANGSD_Total_Sites_down / ANGSD_Total_Sites_orig,
    SNPAD_Total = SNPAD_Total_Sites_down / SNPAD_Total_Sites_orig
  ) %>%
  pivot_longer(
    cols = -c(Individual, Group, MeanCoverage_down),
    names_to = "Method_Metric",
    values_to = "Retention"
  ) %>%
  separate(Method_Metric, into = c("Method", "Metric"), sep = "_")

p_ret_compare <- ggplot(
  ret_df2,
  aes(x = Method, y = Retention, colour = Group)
) +
 geom_boxplot(
    aes(group = interaction(Method, Group)),
    width = 0.7,
    alpha = 0.8,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    colour = "black"
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.7
    ),
    alpha = 0.7
  ) +
    ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 1,
    show.legend = FALSE
  ) +
  facet_grid(Metric ~ Group) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  labs(
    x = "",
    y = "Retention (downsampled / original)",
    title = "Effect of downsampling on total vs heterozygous sites",
    subtitle = "Testing whether heterozygous sites are disproportionately affected"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

print(p_ret_compare)
ggsave("plot_down5x_retention_compare.png", p_ret_compare, width = 8, height = 8)





# QUESTION 4: DOES COVERAGE AFFECT HET SITES DIFFERENTLY FOR THE TWO METHODS?
# RESULTS

# Coverage x Het sites
p_cov_het <- ggplot(
  plot_df,
  aes(x = MeanCoverage, y = Het_Sites, colour = Group)
) +
  geom_point(aes(shape = Group), size = 2.5, alpha = 0.8) +
  geom_smooth(aes(group = Group), method = "lm", se = FALSE, linewidth = 0.8) +
  facet_grid(Treatment~ Method, scales = "free_x") +
  labs(
    x = "Original mean coverage of sample (×)",
    y = "Number of heterozygous sites",
    title = "Effect of downsampling on heterozygous site counts",
    subtitle = "Scaffold 1; Only individuals successfully downsampled to ~5× included. Lines show linear regression fits for each method and treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  ) +
    ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 1,
    show.legend = FALSE
  ) 

print(p_cov_het)
ggsave("plot_coverage_vs_het_sites_by_group_scaf1.png", p_cov_het, width = 9, height = 9)
 

# coverage threshold
p_thresh <- ggplot(
  ret_df,
  aes(x = MeanCoverage_orig, y = Het_retention, colour = Group)
) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +
  facet_grid(Method ~ .) +
   ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 1,
    show.legend = FALSE
  ) +
  labs(
    x = "Original mean coverage of sample (×)",
    y = "Het retention (down / original)",
    title = "Coverage threshold effects on heterozygosity estimation",
    subtitle = "Only individuals successfully downsampled to ~5× included"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 5, linetype = 2, colour = "grey40")

print(p_thresh)
ggsave("plot_coverage_threshold_effects.png", p_thresh, width = 14, height = 8)
 


 # QUESTION 5: Does a larger covreage drop lead to larger het site loss? Does it lead to decreased het estimates?
 # QUESTION: Is ANGSD het increase after downsampling related to original coverage / coverage drop?

coverage_change_df <- df_filt %>%
  transmute(
    Individual,
    Group,
    Year,
    MeanCoverage_orig = as.numeric(MeanCoverage),
    MeanCoverage_down = as.numeric(MeanCoverage_down5x),
    Coverage_drop = MeanCoverage_orig - MeanCoverage_down,
    Coverage_ratio = MeanCoverage_down / MeanCoverage_orig,

    ANGSD_Het_Sites_orig = ANGSD_Het_Sites_orig,
    ANGSD_Het_Sites_down = ANGSD_Het_Sites_down,
    ANGSD_Total_Sites_orig = ANGSD_Total_Sites_orig,
    ANGSD_Total_Sites_down = ANGSD_Total_Sites_down
  ) %>%
  mutate(
    ANGSD_Het_Est_orig = ANGSD_Het_Sites_orig / ANGSD_Total_Sites_orig,
    ANGSD_Het_Est_down = ANGSD_Het_Sites_down / ANGSD_Total_Sites_down,

    Het_sites_change = ANGSD_Het_Sites_down - ANGSD_Het_Sites_orig,
    Het_sites_pct_change = 100 * (ANGSD_Het_Sites_down - ANGSD_Het_Sites_orig) / ANGSD_Het_Sites_orig,

    Het_est_change = ANGSD_Het_Est_down - ANGSD_Het_Est_orig,
    Het_est_pct_change = 100 * (ANGSD_Het_Est_down - ANGSD_Het_Est_orig) / ANGSD_Het_Est_orig,

    Increased_het_est = ANGSD_Het_Est_down > ANGSD_Het_Est_orig
  )

# Change in het estimate
p_angsd_covdrop <- ggplot(
  coverage_change_df,
  aes(x = Coverage_drop, y = Het_est_pct_change, colour = Group)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(shape = Increased_het_est), size = 3, alpha = 0.85) +
  #geom_smooth(aes(group = Group), method = "lm", se = FALSE, linewidth = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 5,
    show.legend = FALSE
  ) +
  labs(
    x = "Coverage drop after downsampling (original coverage - downsampled coverage)",
    y = "Change in ANGSD heterozygosity estimate (%)",
    title = "Is increased ANGSD heterozygosity after downsampling related to coverage drop?",
    subtitle = "Values above zero increased after downsampling to ~5×",
    colour = "Group",
    shape = "ANGSD het estimate increased"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

print(p_angsd_covdrop)

ggsave(
  "plot_angsd_het_est_change_vs_coverage_drop.png",
  p_angsd_covdrop,
  width = 9,
  height = 6,
  dpi = 300
)

# Change in het sites
p_angsd_covdrop_sites <- ggplot(
  coverage_change_df,
  aes(x = Coverage_drop, y = Het_sites_pct_change, colour = Group)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(shape = Increased_het_est), size = 3, alpha = 0.85) +
  #geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = Individual),
    size = 2.7,
    max.overlaps = 3,
    show.legend = FALSE
  ) +
  labs(
    x = "Coverage drop after downsampling (original coverage - downsampled coverage)",
    y = "Change in ANGSD heterozygous sites (%)",
    title = "Is change in ANGSD heterozygous site count related to coverage drop?",
    subtitle = "Values above zero mean more ANGSD het sites after downsampling",
    colour = "Group",
    shape = "ANGSD het estimate increased"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

print(p_angsd_covdrop_sites)

ggsave(
  "plot_angsd_het_sites_change_vs_coverage_drop.png",
  p_angsd_covdrop_sites,
  width = 9,
  height = 6,
  dpi = 300
)

library(patchwork)

p_combined <- p_angsd_covdrop_sites + p_angsd_covdrop + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(p_combined)
ggsave("plot_angsd_cov_drop_combined.png", p_combined, width = 16, height = 8, dpi = 300)







# CLAUDIA - Testing which denominator to use 
# Questions
#1. How different are the two denominators?
#  2. Does denominator choice affect ANGSD and snpAD differently?
#  3. Does denominator choice affect original and downsampled data differently?
#  4. Is the denominator effect associated with coverage?
#  5. Which denominator gives more stable estimates after downsampling?

setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/Analyses/Het/Claudia_denominator_tests")

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# Compare denominator choice for heterozygosity estimates
# Individual callable sites vs Luke same-site denominator
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 0. Make one clean long dataframe
# ============================================================

df_filt2 <- df_filt %>%
  filter(MeanCoverage >= 4)

df_long <- bind_rows(
  df_filt2 %>%
    transmute(
      Individual, Group, Year, MeanCoverage,
      Treatment = "Original",
      Method = "ANGSD",
      Het_sites = ANGSD_Het_Sites_orig,
      Individual_callable_sites = ANGSD_Total_Sites_orig,
      Luke_sites = Luke_Scaf1_Count_orig
    ),
  
  df_filt2 %>%
    transmute(
      Individual, Group, Year, MeanCoverage,
      Treatment = "Original",
      Method = "snpAD",
      Het_sites = SNPAD_Het_Sites_orig,
      Individual_callable_sites = SNPAD_Total_Sites_orig,
      Luke_sites = Luke_Scaf1_Count_orig
    ),
  
  df_filt2 %>%
    transmute(
      Individual, Group, Year, MeanCoverage,
      Treatment = "Downsampled 5x",
      Method = "ANGSD",
      Het_sites = ANGSD_Het_Sites_down,
      Individual_callable_sites = ANGSD_Total_Sites_down,
      Luke_sites = Luke_Scaf1_Count_down
    ),
  
  df_filt2 %>%
    transmute(
      Individual, Group, Year, MeanCoverage,
      Treatment = "Downsampled 5x",
      Method = "snpAD",
      Het_sites = SNPAD_Het_Sites_down,
      Individual_callable_sites = SNPAD_Total_Sites_down,
      Luke_sites = Luke_Scaf1_Count_down
    )
) %>%
  mutate(
    Het_individual_denom = Het_sites / Individual_callable_sites,
    Het_luke_denom = Het_sites / Luke_sites,
    Treatment = factor(Treatment, levels = c("Original", "Downsampled 5x")),
    Method = factor(Method, levels = c("ANGSD", "snpAD"))
  )

# Order individuals by year
sample_order <- df_long %>%
  arrange(Year, Individual) %>%
  pull(Individual) %>%
  unique()

df_long <- df_long %>%
  mutate(
    Individual = factor(Individual, levels = sample_order)
  )



# PLOT 1: DOES DENOM CHANGE HET ESTIMATES?
# One coverage label per individual per facet
cov_labels_p1 <- denom_plot_df %>%
  group_by(Individual, Treatment, Method) %>%
  summarise(
    MeanCoverage = first(MeanCoverage),
    y_pos = max(Het_estimate, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    cov_label = paste0(round(MeanCoverage, 1), "x")
  )

p1 <- ggplot(
  denom_plot_df %>% filter(Treatment == "Original"),
  aes(x = Individual, y = Het_estimate, fill = Denominator)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    data = cov_labels_p1 %>% filter(Treatment == "Original"),
    aes(
      x = Individual,
      y = y_pos,
      label = cov_label,
    ),
    inherit.aes = FALSE,
    angle = 0,
    size = 2.5,
    vjust = 0.5,
    fontface = "bold"
  ) +
  facet_grid(Method ~ ., scales = "fixed") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(
    x = "Individual",
    y = "Heterozygosity estimate",
    fill = "Denominator",
    title = "Does denominator choice change heterozygosity estimates?",
    subtitle = "Comparison of individual callable sites vs Luke same-site denominator"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p1)
ggsave("p1_raw_het_denom_test.png", p1, width = 14, height = 10)

# PLOT 2: Using Luke denominator, do original and downsampled estimates match?
luke_method_df <- df_long %>%
  
  select(
    
    Individual, Group, Year, MeanCoverage,
    
    Treatment, Method,
    
    Het_luke_denom
  )
cov_labels_p2 <- luke_method_df %>%
  group_by(Individual, Method) %>%
  summarise(
    MeanCoverage = first(MeanCoverage),
    y_pos = max(Het_luke_denom, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    cov_label = paste0(round(MeanCoverage, 1), "x")
  )

p2 <- ggplot(
  luke_method_df,
  aes(x = Individual, y = Het_luke_denom, fill = Treatment)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    data = cov_labels_p2,
    aes(
      x = Individual,
      y = y_pos,
      label = cov_label,
      fontface = "bold"
    ),
    inherit.aes = FALSE,
    angle = 0,
    size = 2.5,
    vjust = 0.5
  ) +
  facet_grid(Method ~ ., scales = "fixed") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(
    x = "Individual",
    y = "Heterozygosity estimate",
    fill = "Treatment",
    title = "Original vs downsampled heterozygosity estimates",
    subtitle = "Using Luke same-site denominator only"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
p2
ggsave("p2_down_compare.png", p2, width = 14, height = 10)

# PLOT 3: What is the ANGSD/snpAD ratio using Luke denominator?
ratio_df <- luke_method_df %>%
  pivot_wider(
    names_from = Method,
    values_from = Het_luke_denom
  ) %>%
  mutate(
    ANGSD_snpAD_ratio = ANGSD / snpAD,
    ANGSD_snpAD_difference = ANGSD - snpAD
  )
cov_labels_p3 <- ratio_df %>%
  group_by(Individual) %>%
  summarise(
    MeanCoverage = first(MeanCoverage),
    y_pos = max(ANGSD_snpAD_ratio, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    cov_label = paste0(round(MeanCoverage, 1), "x")
  )

p3 <- ggplot(
  ratio_df,
  aes(x = Individual, y = ANGSD_snpAD_difference, fill = Treatment)
) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    data = cov_labels_p3,
    aes(
      x = Individual,
      y = y_pos,
      label = cov_label
    ),
    inherit.aes = FALSE,
    angle = 0,
    size = 5,
    vjust = 0.5,
    fontface = "bold"
  ) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(
    x = "Individual",
    y = "ANGSD / snpAD heterozygosity difference",
    fill = "Treatment",
    title = "How much heterozygosity does ANGSD capture relative to snpAD?",
    subtitle = "Using Luke same-site denominator only"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
p3
ggsave("p3_ratio.png", p3, width = 14, height = 10)


# PLOT 4: Is the denominator/method/downsampling effect related to coverage?
coverage_denom_df <- df_long %>%
  mutate(
    denom_ratio = Het_luke_denom / Het_individual_denom,
    denom_difference = Het_luke_denom - Het_individual_denom
  )

p4 <- ggplot(
  coverage_denom_df,
  aes(x = MeanCoverage, y = denom_ratio, colour = Method)
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_grid(Method ~ Treatment) +
  theme_bw() +
  labs(
    x = "Original mean coverage",
    y = "Het estimate ratio (Luke denominator / individual denominator)",
    colour = "Method",
    title = "Is denominator effect related to coverage?"
  )
p4
ggsave("p4_cov_denom_effect.png", p4, width = 10, height = 10)

p5 <- ggplot(
  ratio_df,
  aes(x = MeanCoverage, y = ANGSD_snpAD_ratio, colour = Treatment)
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Original mean coverage",
    y = "ANGSD / snpAD heterozygosity ratio",
    colour = "Treatment",
    title = "Is the ANGSD vs snpAD difference related to coverage?",
    subtitle = "Using Luke same-site denominator only"
  )
p5
ggsave("p5_cov_ratio.png", p5, width = 10, height = 10)

# PLOT 6: IS CALLABLE SITE RETENTION DRIVING THIS ESTIMATE?
callable_df <- df_long %>%
  mutate(
    callable_site_ratio = Individual_callable_sites / Luke_sites
  ) %>%
  select(
    Individual,
    Group,
    Year,
    MeanCoverage,
    Treatment,
    Method,
    callable_site_ratio,
    Het_individual_denom,
    Het_luke_denom
  ) %>%
  pivot_longer(
    cols = c(Het_individual_denom, Het_luke_denom),
    names_to = "Denominator",
    values_to = "Het_estimate"
  ) %>%
  mutate(
    Denominator = recode(
      Denominator,
      Het_individual_denom = "Individual callable sites",
      Het_luke_denom = "Luke same-site denominator"
    )
  )

p6 <- ggplot(
  callable_df,
  aes(
    x = MeanCoverage,
    y = Het_estimate,
    colour = Denominator
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(Method ~ Treatment) +
  theme_bw() +
  labs(
    x = "Individual callable sites / Luke sites",
    y = "Heterozygosity estimate",
    colour = "Denominator",
    title = "How does callable-site retention affect heterozygosity estimates?"
  )

p6
ggsave("p6_cov.png", p6, width = 10, height = 10)
