library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(readxl)

# ---- 1. Read in the data ----
# Expects the long-format CSV with columns: Sample, Group, Stage, Het_per_callable, ...
setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP")

df <- read_excel(
  "~/Library/CloudStorage/OneDrive-UniversityofKent/PhD stuff/PhD research/Research projects/01_OBP/OBP_sumstats.xlsx",
  sheet = "vcf_snpad_plotting"
) 


# ---- 2. Set stage order and nicer labels ----
stage_levels <- c(
  "1_Raw_gVCF",
  "2_Filtered_gVCF",
  "3_Strict_miss80",
  "4_Strict_miss100",
  "5_Bisnp_with_singletons",
  "6_Bisnp_no_singletons"
)

stage_labels <- c(
  "Raw gVCF",
  "Filtered gVCF",
  "Strict\n(miss<20%)",
  "Strict\n(miss=0%)",
  "Biallelic SNPs\n(+singletons)",
  "Biallelic SNPs\n(no singletons)"
)

df <- df %>%
  mutate(Stage = factor(Stage, levels = stage_levels, labels = stage_labels))



# ---- 3. Prepare data ----

# Order individuals by year, then by sample name within each year
sample_order <- df %>%
  distinct(Sample, Group, Year) %>%
  arrange(Group, Year, Sample) %>%
  pull(Sample)

df <- df %>%
  mutate(
    Sample = factor(Sample, levels = sample_order)
  )

# ---- 4. Grouped bar plot ----
p <- ggplot(
  #df,
  df %>% filter(Stage %in% c("Biallelic SNPs\n(+singletons)", "Biallelic SNPs\n(no singletons)")),
  #df %>% filter(Stage %in% c("Strict\n(miss<20%)", "Strict\n(miss=0%)")),
  aes(
    x = Sample,
    y = Het_per_callable,
    fill = Stage
  )
) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.8
  ) +
  facet_grid(
    . ~ Group,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    #title = "OBP genome-wide autosomal Heterozygosity across the variant-calling pipeline",
    title = "OBP heterozygosity - Comparing singletons vs no singletons",
    subtitle = "Ordered chronologically. Used the miss20 vcf.",
    x = "Individual",
    y = "Genome-wide heterozygosity",
    fill = "Pipeline stage"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(1.2, "cm"),
    panel.spacing.x = grid::unit(1, "lines")
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

print(p)

# ---- 5. Save the grouped bar plot ----
ggsave(
  filename = "heterozygosity_barplot_by_individual.png",
  plot = p,
  width = 16,
  height = 7,
  dpi = 300
)





# ---- 3. Plot: heterozygosity per individual across stages, faceted by group ----
p <- ggplot(df, aes(x = Stage, y = Het_per_callable, group = Sample, color = Group)) +
  geom_line(alpha = 0.6, linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_color_manual(values = c("Modern" = "#3B6E8F",
                                 "Historical" = "#C9622D")) +
  labs(
    title = "Heterozygosity per callable site across the variant-calling pipeline",
    subtitle = "Each line is one individual, faceted by sample group",
    x = NULL,
    y = "Heterozygous sites / callable sites",
    color = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p)

# ---- 4. Save to file ----
ggsave("heterozygosity_by_stage_faceted.png", p, width = 10, height = 5.5, dpi = 300)
