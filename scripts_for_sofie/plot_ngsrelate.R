# ==============================================================================
# PLOT NgsRelate RESULTS AND LIST RELATED PAIRS
#
# Inputs (produced by run_ngsrelate.sh):
#   <base_name>.ngsrelate
#   <base_name>.ids
#
# Outputs (written to the working directory):
#   <base_name>_relatedness.pdf        heatmap of pairwise theta
#   <base_name>_related_pairs.tsv      all pairs at second-degree or closer
#
# Kinship thresholds (theta), as in Fontsere et al. (2025):
#   first-degree or closer:  theta >= 0.1875
#   second-degree:           0.09375 <= theta < 0.1875
# ==============================================================================

library(Matrix)
library(corrplot)

# ------------------------------------------------------------
# USER SETTINGS (EDIT)
# ------------------------------------------------------------
setwd("/path/to/folder/with/ngsrelate/results")   # EDIT

base_name <- "combined"          # same as NAME in run_ngsrelate.sh
species   <- "little owl"        # used in the plot title

first_degree_min  <- 0.1875
second_degree_min <- 0.09375


# ------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------
rel   <- read.table(paste0(base_name, ".ngsrelate"), header = TRUE)
Names <- scan(paste0(base_name, ".ids"), what = "character", quiet = TRUE)
n     <- length(Names)
N_sites <- rel$nSites[1]

# NgsRelate numbers samples from 0; R numbers from 1.
rel$a <- rel$a + 1L
rel$b <- rel$b + 1L

if (max(c(rel$a, rel$b)) > n) {
  stop("The .ngsrelate file has more samples than the .ids file. ",
       "Check that both come from the same run.")
}

rel$name_a <- Names[rel$a]
rel$name_b <- Names[rel$b]


# ------------------------------------------------------------
# THETA MATRIX + HEATMAP
# ------------------------------------------------------------
# sparseMatrix(symmetric = TRUE) needs every pair stored as row <= column,
# so put the smaller index first.
mat <- sparseMatrix(
  i = pmin(rel$a, rel$b),
  j = pmax(rel$a, rel$b),
  x = rel$theta,
  dims = c(n, n),
  symmetric = TRUE
)
mat <- as.matrix(mat)
diag(mat) <- 0.5                       # an individual's relatedness to itself
dimnames(mat) <- list(Names, Names)

# Shrink labels a little when there are many samples.
label_size <- if (n > 40) 0.5 else 0.95

pdf(paste0(base_name, "_relatedness.pdf"), width = 8, height = 8)
corrplot(
  mat,
  method = "square",
  type   = "lower",
  is.corr = FALSE,
  col.lim = c(0, 0.5),
  tl.col = "black",
  tl.cex = label_size,
  tl.srt = 45
)
mtext(
  paste0("Relatedness among ", species, " individuals (", base_name,
         "; n = ", n, ", sites = ", N_sites, ")"),
  side = 1, line = 4, cex = 1.1
)
dev.off()


# ------------------------------------------------------------
# RELATED PAIRS
# ------------------------------------------------------------
rel$degree <- ifelse(rel$theta >= first_degree_min, "first_degree_or_closer",
              ifelse(rel$theta >= second_degree_min, "second_degree",
                     "unrelated_or_distant"))

related <- rel[rel$degree != "unrelated_or_distant",
               c("name_a", "name_b", "theta", "degree")]
related <- related[order(-related$theta), ]

print(related, row.names = FALSE)

write.table(related, paste0(base_name, "_related_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Individuals appearing in more than one first-degree pair are often the
# best single sample to remove (as with the OBP dataset).
first <- related[related$degree == "first_degree_or_closer", ]
if (nrow(first) > 0) {
  cat("\nNumber of first-degree pairs each individual is in:\n")
  print(sort(table(c(first$name_a, first$name_b)), decreasing = TRUE))
} else {
  cat("\nNo first-degree pairs found.\n")
}
