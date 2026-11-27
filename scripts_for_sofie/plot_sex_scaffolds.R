# ==============================================================================
# PLOT RESULTS FROM sex_scaffolds.sh
#
# Inputs (copy these from OUT_DIR on the cluster into one folder):
#   sex_results.tsv
#   all_normalised_depths.tsv
#   possible_sex_linked_scaffolds.tsv
#   (optional) a sample information file, see SAMPLE_INFO below
#
# Outputs (written to the working directory):
#   1. sex_scaffold_screen.png
#        Every scaffold: mean normalised depth in males (x) vs females (y).
#        Autosomes sit near (1, 1), Z-linked scaffolds near (1, 0.5),
#        W-linked scaffolds near (0, 0.5).
#   2. candidate_scaffold_individual_depths.png
#        Each sample's normalised depth at the candidate scaffolds, by sex.
#   3. normalised_depth_heatmap.png
#        Samples x scaffolds heatmap (control autosome + candidates).
#
# The candidate scaffolds are read from possible_sex_linked_scaffolds.tsv,
# so you do not need to type scaffold names, unless you want to.
# ==============================================================================

library(tidyverse)


# ------------------------------------------------------------
# USER SETTINGS (EDIT)
# ------------------------------------------------------------
setwd("/path/to/folder/with/results")   # EDIT

species <- "Species name"               # EDIT: used in plot titles

# Optional: a tab-separated file with at least the columns
#   Sample   (must match the Sample column in sex_results.tsv)
#   Group    (e.g. historical / modern, or population)
# Points are then coloured by Group. Leave as NA if you do not have one.
SAMPLE_INFO <- NA                        # e.g. "sample_info.tsv"

# Scaffolds to plot. Leave as NULL to use every scaffold classified as
# possible_Z_linked or possible_W_linked (the largest MAX_CANDIDATES).
candidate_scaffolds <- NULL              # e.g. c("scaffold_5", "scaffold_53")
MAX_CANDIDATES <- 8

# Autosomal control scaffold for the heatmap. NULL = the longest scaffold
# classified as autosomal_or_unknown.
control_scaffold <- NULL                 # e.g. "scaffold_1"


# ------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------
sex_results <- read_tsv("sex_results.tsv", show_col_types = FALSE)
all_depths  <- read_tsv("all_normalised_depths.tsv", show_col_types = FALSE)
screen      <- read_tsv("possible_sex_linked_scaffolds.tsv", show_col_types = FALSE,
                        na = c("", "NA"))

sex_levels <- c("male_ZZ", "female_ZW", "ambiguous", "low_coverage")

sex_info <- sex_results %>%
  select(Sample, Predicted_sex) %>%
  mutate(Predicted_sex = factor(Predicted_sex, levels = sex_levels))

# Add the optional sample groups (e.g. historical / modern).
if (!is.na(SAMPLE_INFO)) {
  groups <- read_tsv(SAMPLE_INFO, show_col_types = FALSE) %>%
    select(Sample, Group)
  sex_info <- left_join(sex_info, groups, by = "Sample")
  missing <- sex_info$Sample[is.na(sex_info$Group)]
  if (length(missing) > 0) {
    warning("No Group in SAMPLE_INFO for: ", paste(missing, collapse = ", "))
  }
} else {
  sex_info$Group <- "all samples"
}

plot_data <- all_depths %>% left_join(sex_info, by = "Sample")

cat("Samples per predicted sex:\n")
print(table(sex_info$Predicted_sex))


# ------------------------------------------------------------
# CHOOSE SCAFFOLDS TO PLOT
# ------------------------------------------------------------
if (is.null(candidate_scaffolds)) {
  candidate_scaffolds <- screen %>%
    filter(Classification %in% c("possible_Z_linked", "possible_W_linked")) %>%
    arrange(desc(Length)) %>%
    slice_head(n = MAX_CANDIDATES) %>%
    pull(Scaffold)
}
if (length(candidate_scaffolds) == 0) {
  stop("No candidate scaffolds found in possible_sex_linked_scaffolds.tsv. ",
       "Set candidate_scaffolds by hand.")
}

if (is.null(control_scaffold)) {
  control_scaffold <- screen %>%
    filter(Classification == "autosomal_or_unknown") %>%
    arrange(desc(Length)) %>%
    slice_head(n = 1) %>%
    pull(Scaffold)
}

cat("\nCandidate scaffolds:", paste(candidate_scaffolds, collapse = ", "), "\n")
cat("Control scaffold:   ", control_scaffold, "\n")


# ------------------------------------------------------------
# PLOT 1: GENOME-WIDE SCREEN (males vs females, every scaffold)
# ------------------------------------------------------------
screen_plot <- ggplot(screen, aes(x = Male_mean, y = Female_mean,
                                  colour = Classification)) +
  geom_point(aes(size = Length), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +     # autosomal
  geom_abline(slope = 0.5, intercept = 0, linetype = "dotted") +   # Z-linked
  scale_size_continuous(range = c(0.5, 5), guide = "none") +
  scale_colour_manual(values = c(autosomal_or_unknown = "grey60",
                                 possible_Z_linked = "#D55E00",
                                 possible_W_linked = "#0072B2")) +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1.5)) +
  labs(
    x = "Mean normalised depth in males",
    y = "Mean normalised depth in females",
    title = paste("Sex-linked scaffold screen:", species),
    subtitle = "Dashed: autosomal (female = male). Dotted: Z-linked (female = 0.5 x male). Point size = scaffold length."
  ) +
  theme_bw()

ggsave("sex_scaffold_screen.png", screen_plot, width = 8, height = 7, dpi = 300)


# ------------------------------------------------------------
# PLOT 2: INDIVIDUAL DEPTHS AT CANDIDATE SCAFFOLDS
# ------------------------------------------------------------
candidate_depths <- plot_data %>%
  filter(Scaffold %in% candidate_scaffolds,
         Predicted_sex %in% c("male_ZZ", "female_ZW")) %>%
  mutate(Scaffold = factor(Scaffold, levels = candidate_scaffolds))

candidate_plot <- ggplot(candidate_depths,
                         aes(x = Predicted_sex, y = Normalised_depth)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(colour = Group), width = 0.12, height = 0,
              size = 2, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  facet_wrap(~ Scaffold, nrow = 1) +
  labs(
    x = "Predicted sex",
    y = "Normalised scaffold depth",
    colour = NULL,
    title = paste("Individual depth at candidate sex-linked scaffolds:", species),
    subtitle = "Z-linked: males ~1, females ~0.5. W-linked: males ~0, females ~0.5."
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("candidate_scaffold_individual_depths.png", candidate_plot,
       width = max(6, 2.2 * length(candidate_scaffolds)), height = 5, dpi = 300)


# ------------------------------------------------------------
# PLOT 3: HEATMAP (all samples, including ambiguous / low_coverage)
# ------------------------------------------------------------
heatmap_scaffolds <- c(control_scaffold, candidate_scaffolds)

heatmap_data <- plot_data %>%
  filter(Scaffold %in% heatmap_scaffolds) %>%
  mutate(
    Scaffold = factor(Scaffold, levels = heatmap_scaffolds),
    Sample   = factor(Sample, levels = rev(sex_results$Sample))
  )

depth_heatmap <- ggplot(heatmap_data,
                        aes(x = Scaffold, y = Sample, fill = Normalised_depth)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1.5),
                       oob = scales::squish, breaks = c(0, 0.5, 1, 1.5),
                       name = "Normalised\ndepth") +
  facet_grid(Predicted_sex ~ ., scales = "free_y", space = "free_y",
             switch = "y", drop = TRUE) +
  labs(
    x = "Scaffold", y = "Sample",
    title = paste("Normalised scaffold depth across", species, "samples"),
    subtitle = paste0("Control autosomal scaffold (", control_scaffold,
                      ") and candidate sex-linked scaffolds")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = if (nrow(sex_results) > 50) 6 else 10)
  )

ggsave("normalised_depth_heatmap.png", depth_heatmap,
       width = max(6, 1.2 * length(heatmap_scaffolds) + 4),
       height = max(5, 0.2 * nrow(sex_results) + 2), dpi = 300)

cat("\nSaved: sex_scaffold_screen.png, candidate_scaffold_individual_depths.png,",
    "normalised_depth_heatmap.png\n")
