# ==============================================================================
# DETECT TECHNICAL OUTLIER SAMPLES FROM A SPREADSHEET OF QUALITY METRICS
#
# What it does
#   1. Reads a spreadsheet: one row per sample, one column per quality metric
#      (e.g. coverage, damage, fragment length, duplication rate).
#   2. Checks the data (missing values, constant columns, skewed columns).
#   3. Reports which metrics repeat each other (Spearman rho, VIF).
#   4. Within each group (e.g. historical, modern), asks two questions:
#        a. Is the sample extreme on any ONE metric?        (robust z-score)
#        b. Is the sample's COMBINATION of metrics unusual? (leave-one-out
#           Mahalanobis distance, with an F test)
#   5. Draws a PCA of each group so you can see where the flagged samples sit
#      and which metrics push them there.
#   6. Writes one summary table: every sample, whether it was flagged, and why.
#
# Flagged samples are candidates to INSPECT, not automatic exclusions.
#
# Input spreadsheet (see quality_metrics_TEMPLATE.xlsx)
#   Sheet "metrics":   one column "sample", optional column "group", then one
#                      numeric column per metric.
#   Sheet "variables": (optional) one row per metric with columns
#                        variable   column name in the metrics sheet
#                        include    yes / no
#                        direction  higher_is_better / lower_is_better / unknown
#                        transform  none / log
#                      Without this sheet, every numeric column is used,
#                      untransformed, with direction "unknown".
#   A .csv file also works (metrics only, no variables sheet).
#
# Packages: readxl, ggplot2
# ==============================================================================

library(readxl)
library(ggplot2)


# ============================================================
# USER SETTINGS (EDIT)
# ============================================================
setwd("/path/to/project/folder")                 # EDIT

DATA_FILE       <- "quality_metrics.xlsx"         # .xlsx or .csv
METRICS_SHEET   <- "metrics"
VARIABLES_SHEET <- "variables"                    # ignored for .csv
OUT_DIR         <- "outlier_output"

SAMPLE_COL <- "sample"
GROUP_COL  <- "group"   # set to NA if all samples belong to one group

# TRUE  = compare each sample with the other samples in its OWN group
#         (recommended for temporal data: a historical sample is only
#         unusual if it is unusual for a historical sample).
# FALSE = one comparison across all samples together.
WITHIN_GROUP <- TRUE

# Outlier thresholds
P_CUTOFF <- 0.01   # combined (Mahalanobis) outlier if p < this
Z_CUTOFF <- 3.5    # single-metric outlier if |robust z| > this (3.5 is a
                   # common cut-off for median/MAD-based z-scores)

# The combined test is "leave one out": each sample is compared with the
# centre and spread of the OTHER samples in its group, so a bad sample
# cannot hide itself by pulling the group average towards it. The p-value
# comes from an F test, which is exact for small groups (if the metrics are
# roughly normal within the group; use transform = log for skewed ones).

# Redundancy checks
MAX_RHO  <- 0.7    # pairs at or above this |rho| are reported
WARN_VIF <- 5      # VIF above this is reported
# Drop metrics above this VIF from the combined test, worst first.
# Default NA = report only. The Mahalanobis distance already allows for
# correlated metrics, and with few samples one or two outliers can create
# a high VIF on their own, so dropping can remove the very metric that
# reveals them. Set e.g. DROP_VIF <- 30 only if the run fails because two
# metrics are almost identical.
DROP_VIF <- NA

# Smallest group size for the combined test, per metric used.
# With 4 metrics and MIN_N_PER_VAR = 3, a group needs at least 12 samples.
MIN_N_PER_VAR <- 3


# ============================================================
# 0. HELPERS
# ============================================================
dir.create(OUT_DIR, showWarnings = FALSE)
log_lines <- c(paste("Run:", Sys.time()), paste("Data:", DATA_FILE))
note <- function(...) {
  msg <- paste0(...)
  message(msg)
  log_lines <<- c(log_lines, msg)
}

# Robust z-score: (value - median) / MAD. Falls back to mean/SD if the MAD
# is zero (e.g. when more than half the samples share one value).
robust_z <- function(x) {
  m <- median(x, na.rm = TRUE)
  s <- mad(x, na.rm = TRUE)
  if (is.na(s) || s == 0) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE) }
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

# VIF for each column: 1 / (1 - R^2), where R^2 comes from predicting the
# column from all the other columns.
vif <- function(d) {
  if (ncol(d) < 2) return(setNames(1, names(d)))
  sapply(names(d), function(v) {
    r2 <- summary(lm(reformulate(setdiff(names(d), v), response = v), data = d))$r.squared
    if (r2 >= 1) Inf else 1 / (1 - r2)
  })
}


# ============================================================
# 1. LOAD AND CHECK THE DATA
# ============================================================
is_xlsx <- grepl("\\.xlsx?$", DATA_FILE, ignore.case = TRUE)
dat <- if (is_xlsx) as.data.frame(read_excel(DATA_FILE, sheet = METRICS_SHEET)) else
                    read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
names(dat) <- trimws(names(dat))

if (!SAMPLE_COL %in% names(dat)) stop("No column called '", SAMPLE_COL, "' in the metrics sheet")
dat[[SAMPLE_COL]] <- trimws(as.character(dat[[SAMPLE_COL]]))
if (anyDuplicated(dat[[SAMPLE_COL]])) {
  stop("Duplicate sample names: ",
       paste(unique(dat[[SAMPLE_COL]][duplicated(dat[[SAMPLE_COL]])]), collapse = ", "))
}

if (is.na(GROUP_COL) || !GROUP_COL %in% names(dat)) {
  dat$.group <- "all"
} else {
  dat$.group <- trimws(as.character(dat[[GROUP_COL]]))
  if (any(is.na(dat$.group) | dat$.group == "")) stop("Some samples have no group")
}
if (!WITHIN_GROUP) dat$.group_fit <- "all" else dat$.group_fit <- dat$.group

# Variable settings: from the "variables" sheet if there is one.
has_var_sheet <- is_xlsx && VARIABLES_SHEET %in% excel_sheets(DATA_FILE)
if (has_var_sheet) {
  vars_tab <- as.data.frame(read_excel(DATA_FILE, sheet = VARIABLES_SHEET))
  names(vars_tab) <- tolower(trimws(names(vars_tab)))
  vars_tab$variable <- trimws(vars_tab$variable)
  unknown <- setdiff(vars_tab$variable, names(dat))
  if (length(unknown)) stop("In the variables sheet but not the metrics sheet: ",
                            paste(unknown, collapse = ", "))
  vars_tab <- vars_tab[tolower(vars_tab$include) %in% c("yes", "y", "true", "1"), ]
} else {
  cand <- setdiff(names(dat), c(SAMPLE_COL, GROUP_COL, ".group", ".group_fit"))
  vars_tab <- data.frame(variable = cand, direction = "unknown", transform = "none")
}
vars_tab$direction[is.na(vars_tab$direction)] <- "unknown"
vars_tab$transform[is.na(vars_tab$transform)] <- "none"
VARS <- vars_tab$variable

# Make sure every metric is numeric.
for (v in VARS) {
  x <- suppressWarnings(as.numeric(dat[[v]]))
  bad <- !is.na(dat[[v]]) & is.na(x)
  if (any(bad)) note("WARNING: non-numeric values in ", v, " set to NA for: ",
                     paste(dat[[SAMPLE_COL]][bad], collapse = ", "))
  dat[[v]] <- x
}

# Log-transform metrics marked "log" (e.g. coverage, which is usually skewed).
for (v in vars_tab$variable[vars_tab$transform == "log"]) {
  if (any(dat[[v]] <= 0, na.rm = TRUE)) {
    note("WARNING: ", v, " has values <= 0; using log10(x + 1) instead of log10(x)")
    dat[[v]] <- log10(dat[[v]] + 1)
  } else dat[[v]] <- log10(dat[[v]])
}

note(nrow(dat), " samples; ", length(VARS), " metrics: ", paste(VARS, collapse = ", "))
note("Samples per group: ", paste(names(table(dat$.group)), table(dat$.group),
                                  sep = " = ", collapse = "; "))
miss <- colSums(is.na(dat[VARS]))
if (any(miss > 0)) note("Missing values: ", paste(names(miss)[miss > 0], miss[miss > 0],
                                                  sep = " = ", collapse = "; "))

# Flag strongly skewed metrics that are not log-transformed.
for (v in vars_tab$variable[vars_tab$transform != "log"]) {
  x <- dat[[v]][!is.na(dat[[v]])]
  if (length(x) > 5 && all(x > 0) && max(x) / median(x) > 5) {
    note("NOTE: ", v, " is strongly skewed (max is > 5 x median). ",
         "Consider transform = log in the variables sheet.")
  }
}


# ============================================================
# 2. ANALYSE EACH GROUP
# ============================================================
results   <- list()
redund    <- list()
pca_plots <- list()

for (g in unique(dat$.group_fit)) {
  note("\n== Group: ", g, " ==")
  dg <- dat[dat$.group_fit == g, ]

  # --- Metrics usable in this group: not constant, not all missing
  use <- VARS[sapply(VARS, function(v) {
    x <- dg[[v]]; sum(!is.na(x)) >= 3 && sd(x, na.rm = TRUE) > 0
  })]
  dropped_const <- setdiff(VARS, use)
  if (length(dropped_const)) note("Not usable in this group (constant, missing, or < 3 values): ",
                                  paste(dropped_const, collapse = ", "))
  if (length(use) == 0) {
    note("Group skipped: no usable metrics (", nrow(dg), " samples). ",
         "Each group needs at least 3 samples with varying values.")
    next
  }

  # --- Single-metric check: robust z for every sample and metric
  Z <- sapply(use, function(v) robust_z(dg[[v]]))
  Z <- matrix(Z, nrow = nrow(dg), dimnames = list(dg[[SAMPLE_COL]], use))

  # --- Redundancy: correlations and VIF, on complete samples
  cc <- complete.cases(dg[, use, drop = FALSE])
  d  <- dg[cc, use, drop = FALSE]
  rho <- suppressWarnings(cor(d, method = "spearman"))
  hi <- which(upper.tri(rho) & abs(rho) >= MAX_RHO, arr.ind = TRUE)
  if (nrow(hi)) note("Pairs with |rho| >= ", MAX_RHO, ": ",
                     paste0(rownames(rho)[hi[, 1]], "-", colnames(rho)[hi[, 2]],
                            " (", round(rho[hi], 2), ")", collapse = "; "))

  multi <- use
  dropped_vif <- character(0)
  repeat {
    if (is.na(DROP_VIF) || length(multi) < 2 || nrow(d) <= length(multi) + 1) break
    vv <- vif(d[, multi, drop = FALSE])
    if (max(vv) <= DROP_VIF) break
    worst <- names(which.max(vv))
    note("Dropped ", worst, " from the combined test (VIF = ", round(max(vv), 1), " > ", DROP_VIF, ")")
    dropped_vif <- c(dropped_vif, worst)
    multi <- setdiff(multi, worst)
  }
  vv <- if (length(multi) >= 2 && nrow(d) > length(multi) + 1) vif(d[, multi, drop = FALSE]) else
        setNames(rep(NA, length(multi)), multi)
  if (any(vv > WARN_VIF, na.rm = TRUE)) note("VIF above ", WARN_VIF, ": ",
      paste(names(vv)[vv > WARN_VIF], round(vv[vv > WARN_VIF], 1), sep = " = ", collapse = "; "))

  redund[[g]] <- data.frame(group = g, variable = use,
                            VIF = round(vv[use], 2),
                            used_in_combined_test = use %in% multi,
                            row.names = NULL)

  # Correlation heatmap for this group
  rd <- as.data.frame(as.table(rho)); names(rd) <- c("v1", "v2", "rho")
  p_rho <- ggplot(rd, aes(v1, v2, fill = rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
    labs(title = paste("Spearman correlations:", g),
         subtitle = paste0("|rho| >= ", MAX_RHO, " = the two metrics carry similar information"),
         x = NULL, y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(OUT_DIR, paste0("correlations_", g, ".png")), p_rho,
         width = 6.5, height = 5.5, dpi = 150)

  # --- Combined check: Mahalanobis distance on the remaining metrics
  k  <- length(multi)
  dm <- dg[cc, multi, drop = FALSE]
  d2 <- rep(NA_real_, nrow(dg)); pval <- rep(NA_real_, nrow(dg))
  enough <- k >= 2 && nrow(dm) >= max(MIN_N_PER_VAR * k, k + 3)
  if (!enough) {
    note("Combined test skipped: ", nrow(dm), " complete samples for ", k,
         " metrics (needs ", max(MIN_N_PER_VAR * k, k + 3),
         "). Only single-metric checks are used.")
  } else {
    X <- as.matrix(dm)
    n <- nrow(X); m <- n - 1
    # Leave one out: distance of sample i from the mean and covariance of
    # the other m = n - 1 samples. For a sample from the same population,
    #   (m / (m + 1)) * d2  ~  (k (m - 1) / (m - k)) * F(k, m - k)
    loo <- sapply(seq_len(n), function(i) {
      rest <- X[-i, , drop = FALSE]
      tryCatch(mahalanobis(X[i, ], colMeans(rest), cov(rest)), error = function(e) NA_real_)
    })
    if (anyNA(loo)) note("WARNING: covariance could not be inverted for some samples. ",
                         "Two metrics are probably near-identical; set DROP_VIF (e.g. 30).")
    Fstat <- (m / (m + 1)) * loo * (m - k) / (k * (m - 1))
    d2[cc]   <- loo
    pval[cc] <- pf(Fstat, k, m - k, lower.tail = FALSE)

    # PCA for display only (the distance above does not depend on it)
    p   <- prcomp(dm, center = TRUE, scale. = TRUE)
    pct <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
    sc  <- data.frame(sample = dg[[SAMPLE_COL]][cc], PC1 = p$x[, 1], PC2 = p$x[, 2],
                      outlier = pval[cc] < P_CUTOFF)
    kk  <- 0.8 * max(abs(sc[, c("PC1", "PC2")])) / max(abs(p$rotation[, 1:2]))
    arr <- data.frame(label = rownames(p$rotation),
                      x = p$rotation[, 1] * kk, y = p$rotation[, 2] * kk)
    pca_plots[[g]] <- ggplot(sc, aes(PC1, PC2)) +
      geom_segment(data = arr, aes(x = 0, y = 0, xend = x, yend = y), inherit.aes = FALSE,
                   colour = "#0072B2", arrow = arrow(length = unit(0.2, "cm"))) +
      geom_text(data = arr, aes(x * 1.15, y * 1.15, label = label), inherit.aes = FALSE,
                colour = "#0072B2", size = 3.2, fontface = "bold") +
      geom_point(aes(colour = outlier), size = 3) +
      geom_text(aes(label = sample, colour = outlier), vjust = -1, size = 2.8,
                check_overlap = TRUE, show.legend = FALSE) +
      scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "red"),
                          labels = c(`FALSE` = "not flagged", `TRUE` = "combined outlier")) +
      labs(title = paste("PCA of quality metrics:", g),
           subtitle = paste0("PC1 ", pct[1], "%, PC2 ", pct[2], "% | n = ", nrow(sc)),
           caption = "Arrows: metrics. Same direction = redundant; right angle = independent.",
           x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"),
           colour = NULL) +
      theme_bw()
    ggsave(file.path(OUT_DIR, paste0("pca_", g, ".png")), pca_plots[[g]],
           width = 7.5, height = 5.5, dpi = 150)
  }

  # --- Summarise each sample: which metrics are extreme, and in which direction
  dir_of <- setNames(vars_tab$direction, vars_tab$variable)
  describe <- function(v, z) {
    side <- if (z > 0) "high" else "low"
    q <- switch(dir_of[[v]],
                higher_is_better = if (z > 0) "better than group" else "worse than group",
                lower_is_better  = if (z < 0) "better than group" else "worse than group",
                "direction unknown")
    sprintf("%s %s (z = %+.1f, %s)", v, side, z, q)
  }
  for (i in seq_len(nrow(dg))) {
    z <- Z[i, ]
    flagged <- names(z)[abs(z) > Z_CUTOFF]
    top <- names(which.max(abs(z)))
    results[[length(results) + 1]] <- data.frame(
      sample = dg[[SAMPLE_COL]][i],
      group  = dg$.group[i],
      compared_with = g,
      metrics_in_combined_test = if (enough) paste(multi, collapse = ", ") else "",
      mahalanobis_d2 = round(d2[i], 2),
      p_value = signif(pval[i], 3),
      combined_outlier = !is.na(pval[i]) && pval[i] < P_CUTOFF,
      n_single_flags = length(flagged),
      single_metric_flags = if (length(flagged)) paste(sapply(flagged, function(v) describe(v, z[[v]])),
                                                       collapse = "; ") else "",
      most_extreme_metric = describe(top, z[[top]]),
      stringsAsFactors = FALSE)
  }
}


# ============================================================
# 3. SUMMARY TABLE AND OVERVIEW PLOTS
# ============================================================
if (length(results) == 0) stop("No group could be analysed. See the messages above.")
res <- do.call(rbind, results)
res$verdict <- ifelse(res$combined_outlier & res$n_single_flags > 0, "INSPECT: both checks",
               ifelse(res$combined_outlier, "inspect: unusual combination",
               ifelse(res$n_single_flags > 0, "inspect: extreme on one metric", "ok")))
res <- res[order(res$verdict == "ok", -res$mahalanobis_d2), ]
write.csv(res, file.path(OUT_DIR, "outlier_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, redund), file.path(OUT_DIR, "redundancy_checks.csv"), row.names = FALSE)

# Each metric on its own: one panel per metric, one point per sample
long <- do.call(rbind, lapply(VARS, function(v)
  data.frame(sample = dat[[SAMPLE_COL]], group = dat$.group, metric = v, value = dat[[v]])))
long$verdict <- res$verdict[match(long$sample, res$sample)]
long$flagged <- long$verdict != "ok"
p_metrics <- ggplot(long, aes(group, value)) +
  geom_boxplot(outlier.shape = NA, colour = "grey60", width = 0.5) +
  geom_jitter(aes(colour = flagged), width = 0.1, height = 0, size = 2) +
  geom_text(data = subset(long, flagged), aes(label = sample), colour = "red",
            size = 2.5, hjust = -0.15, check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "red"),
                      labels = c(`FALSE` = "ok", `TRUE` = "flagged")) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "Each quality metric, by group",
       subtitle = "Red = flagged by at least one check (see outlier_summary.csv)",
       x = NULL, y = "Value (log10 where transform = log)", colour = NULL) +
  theme_bw()
ggsave(file.path(OUT_DIR, "metrics_by_group.png"), p_metrics,
       width = 10, height = 7, dpi = 150)

# Ranked combined distance, with the cut-off for each group
rk <- res[!is.na(res$mahalanobis_d2), ]
if (nrow(rk)) {
  rk$k <- lengths(strsplit(rk$metrics_in_combined_test, ", "))
  rk$m <- ave(rk$mahalanobis_d2, rk$compared_with, FUN = length) - 1
  # Distance at which p = P_CUTOFF (the F test above, turned back into a distance)
  rk$cut <- sqrt(qf(1 - P_CUTOFF, rk$k, rk$m - rk$k) *
                 rk$k * (rk$m - 1) / (rk$m - rk$k) * (rk$m + 1) / rk$m)
  rk$dist <- sqrt(rk$mahalanobis_d2)
  rk$sample <- factor(rk$sample, levels = rk$sample[order(rk$dist)])
  p_rank <- ggplot(rk, aes(dist, sample, fill = combined_outlier)) +
    geom_col() +
    geom_vline(aes(xintercept = cut), linetype = "dashed", colour = "red") +
    facet_wrap(~ compared_with, scales = "free_y") +
    scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "red"),
                      labels = c(`FALSE` = "not flagged", `TRUE` = "combined outlier")) +
    labs(title = "How unusual is each sample's combination of metrics?",
         subtitle = paste0("Distance from the centre of its group. Dashed = p < ", P_CUTOFF, " cut-off"),
         x = "Mahalanobis distance", y = NULL, fill = NULL) +
    theme_bw()
  ggsave(file.path(OUT_DIR, "distance_ranking.png"), p_rank,
         width = 9, height = max(4, 0.22 * nrow(rk) + 2), dpi = 150)
}

note("\n== Result ==")
note(sum(res$verdict != "ok"), " of ", nrow(res), " samples flagged for inspection:")
for (i in which(res$verdict != "ok")) note("  ", res$sample[i], " (", res$group[i], "): ",
                                          res$verdict[i], " | ", res$most_extreme_metric[i])

log_lines <- c(log_lines, "", "Settings:",
               paste0("  WITHIN_GROUP = ", WITHIN_GROUP,
                      ", P_CUTOFF = ", P_CUTOFF, ", Z_CUTOFF = ", Z_CUTOFF,
                      ", DROP_VIF = ", DROP_VIF, ", MIN_N_PER_VAR = ", MIN_N_PER_VAR))
writeLines(log_lines, file.path(OUT_DIR, "run_log.txt"))
cat("\nOutputs in", normalizePath(OUT_DIR), "\n")
