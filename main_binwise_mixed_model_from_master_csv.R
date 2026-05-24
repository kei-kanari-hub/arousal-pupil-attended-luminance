# ============================================================
# Main bin-wise mixed-effects analysis from the filtered master CSV
#
# This script is designed for repository use. It requires only:
#   all_subjects_master_table_main_analysis_filtered.csv
#
# Usage:
#   Rscript main_binwise_mixed_model_from_master_csv.R \
#     data/all_subjects_master_table_main_analysis_filtered.csv \
#     analysis_outputs/main_binwise_mixed_model
#
# If no arguments are supplied, the script searches for the input CSV in:
#   1) data/all_subjects_master_table_main_analysis_filtered.csv
#   2) all_subjects_master_table_main_analysis_filtered.csv
# and writes output to:
#   analysis_outputs/main_binwise_mixed_model/run_YYYYMMDD_HHMMSS
#
# Outputs:
#   - figure_main_pupil_timecourse.pdf
#   - figure_main_mixed_beta_timecourses_4panel.pdf
#   - figure_main_mixed_beta_timecourses_1_arousal.pdf
#   - figure_main_mixed_beta_timecourses_2_valence.pdf
#   - figure_main_mixed_beta_timecourses_3_block.pdf
#   - figure_main_mixed_beta_timecourses_4_arousal_block.pdf
#   - figure_main_mixed_beta_timecourses_5_valence_block.pdf
#   - grand_mean_pupil_main.csv
#   - binwise_mixed_model_results_raw.csv
#   - binwise_mixed_model_results_fdr.csv
#   - trial_level_affect_zscores.csv
#   - trial_level_affect_zscore_summary.csv
#   - trial_level_affect_zscore_by_subject_summary.csv
#   - trial_time_bin_counts.csv
#   - trial_time_bin_count_summary.csv
#   - trial_time_bin_count_distribution.csv
#   - trial_time_bin_count_problem_trials.csv
#   - run_info.txt
#
# Notes:
#   - This script does not recreate supplementary long-window pupil or OKN
#     visualizations because those require separate full time-series or OKN
#     files. It reproduces the main bin-wise pupil analysis from the filtered
#     master table only.
#   - Arousal and valence are z-standardized within participant after reducing
#     the data to one row per subject x trial.
#   - Block is effect-coded: block 1 = -0.5, block 2 = +0.5.
#   - P values are obtained from lmerTest and FDR-corrected across time bins
#     separately for each effect using the Benjamini-Hochberg procedure.
# ============================================================

rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lmerTest)
  library(grid)
})

# ============================================================
# 1. Command-line arguments and paths
# ============================================================
args <- commandArgs(trailingOnly = TRUE)

find_default_master_file <- function() {
  candidates <- c(
    file.path("data", "all_subjects_master_table_main_analysis_filtered.csv"),
    "all_subjects_master_table_main_analysis_filtered.csv"
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  normalizePath(hit[1], mustWork = TRUE)
}

MASTER_FILE <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  find_default_master_file()
}

if (is.na(MASTER_FILE) || !file.exists(MASTER_FILE)) {
  stop(
    paste0(
      "Input CSV not found. Provide the path as the first argument, e.g.\n",
      "  Rscript 05_main_binwise_mixed_model_from_master_csv.R ",
      "data/all_subjects_master_table_main_analysis_filtered.csv"
    ),
    call. = FALSE
  )
}
MASTER_FILE <- normalizePath(MASTER_FILE, mustWork = TRUE)

OUTPUT_ROOT <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("analysis_outputs", "main_binwise_mixed_model")
}

RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("run_", RUN_TAG))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUTPUT_DIR <- normalizePath(OUTPUT_DIR, mustWork = TRUE)

# ============================================================
# 2. Settings
# ============================================================
ALPHA_Q <- 0.05
DEPENDENT_VAR <- "pupil_pct_change_sound_pre500to0_bin"
MIN_ROWS_PER_BIN <- 10
MIN_SUBJECTS_PER_BIN <- 2

MODEL_FORMULA <- as.formula(
  paste0(
    DEPENDENT_VAR,
    " ~ arousal_z + valence_z + block_ec + ",
    "arousal_z:block_ec + valence_z:block_ec + (1 | subject)"
  )
)

# ============================================================
# 3. Helper functions
# ============================================================
safe_zscore <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mu) / sdv
}

ci95_t <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= 1) {
    return(NA_real_)
  }
  sem <- sd(x) / sqrt(n)
  qt(0.975, df = n - 1) * sem
}

get_coef_stats <- function(coef_df, term_name) {
  out <- list(
    estimate = NA_real_,
    se = NA_real_,
    df = NA_real_,
    t = NA_real_,
    p = NA_real_
  )
  if (is.null(coef_df) || nrow(coef_df) == 0) {
    return(out)
  }
  rn <- rownames(coef_df)
  idx <- which(rn == term_name)
  if (length(idx) == 0) {
    return(out)
  }
  idx <- idx[1]
  out$estimate <- suppressWarnings(as.numeric(coef_df[idx, "Estimate"]))
  out$se <- suppressWarnings(as.numeric(coef_df[idx, "Std. Error"]))
  out$df <- suppressWarnings(as.numeric(coef_df[idx, "df"]))
  out$t <- suppressWarnings(as.numeric(coef_df[idx, "t value"]))
  out$p <- suppressWarnings(as.numeric(coef_df[idx, "Pr(>|t|)"]))
  out
}

interpret_direction <- function(beta) {
  if (!is.finite(beta)) return("NA")
  if (beta > 0) return("positive")
  if (beta < 0) return("negative")
  "zero"
}

summarize_effect_pattern <- function(sig_arousal, sig_valence, sig_block,
                                     sig_arousal_block, sig_valence_block) {
  tags <- character(0)
  if (isTRUE(sig_arousal)) tags <- c(tags, "arousal_main")
  if (isTRUE(sig_valence)) tags <- c(tags, "valence_main")
  if (isTRUE(sig_block)) tags <- c(tags, "block_main")
  if (isTRUE(sig_arousal_block)) tags <- c(tags, "arousal_x_block")
  if (isTRUE(sig_valence_block)) tags <- c(tags, "valence_x_block")
  if (length(tags) == 0) return("none_q<.05")
  paste(tags, collapse = ";")
}

plot_two_ci_lines <- function(df, xcol,
                              y1, ci1, y2, ci2,
                              label1 = "Block 1: black-to-white",
                              label2 = "Block 2: white-to-black",
                              title_txt,
                              ylab_txt,
                              xlab_txt = "Time from sound onset (s)") {
  dd <- bind_rows(
    df %>% transmute(x = .data[[xcol]], mean = .data[[y1]], ci = .data[[ci1]], block = label1),
    df %>% transmute(x = .data[[xcol]], mean = .data[[y2]], ci = .data[[ci2]], block = label2)
  )

  ggplot(dd, aes(x = x, y = mean, color = block, fill = block, linetype = block)) +
    geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(x = xlab_txt, y = ylab_txt, title = title_txt, color = NULL, fill = NULL, linetype = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

plot_beta_panel <- function(df, beta_col, q_col, title_txt, ylab_txt,
                            alpha_q = 0.05, show_title = TRUE) {
  dd <- df %>%
    transmute(
      time_bin_s = time_bin_s,
      beta = .data[[beta_col]],
      q = .data[[q_col]]
    )

  dd2 <- dd[is.finite(dd$time_bin_s) & is.finite(dd$beta), , drop = FALSE]
  sig_df <- dd2 %>% filter(is.finite(q) & q < alpha_q)
  rug_y <- if (nrow(dd2) > 0) min(dd2$beta, na.rm = TRUE) else 0

  p <- ggplot(dd2, aes(x = time_bin_s, y = beta)) +
    geom_line(color = "black", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(
      data = sig_df,
      aes(x = time_bin_s, y = rug_y),
      inherit.aes = FALSE,
      size = 1.4
    ) +
    labs(x = "Time from sound onset (s)", y = ylab_txt) +
    theme_bw(base_size = 11)

  if (show_title) {
    p <- p + labs(title = title_txt)
  }
  p
}

save_four_panel_pdf <- function(p1, p2, p3, p4, filename, width = 8.5, height = 7.2) {
  pdf(filename, width = width, height = height)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  dev.off()
}

# ============================================================
# 4. Load filtered master table
# ============================================================
M <- read_csv(MASTER_FILE, show_col_types = FALSE)

reqM <- c(
  "subject", "trial", "block", "arousal", "valence",
  "time_bin_s", DEPENDENT_VAR
)
missingM <- setdiff(reqM, names(M))
if (length(missingM) > 0) {
  stop(sprintf("Master table missing columns: %s", paste(missingM, collapse = ", ")))
}

M <- M %>%
  mutate(
    subject = as.character(subject),
    trial = as.numeric(trial),
    block = as.numeric(block),
    arousal = as.numeric(arousal),
    valence = as.numeric(valence),
    time_bin_s = as.numeric(time_bin_s),
    pupil_pct_change_sound_pre500to0_bin = as.numeric(.data[[DEPENDENT_VAR]])
  ) %>%
  filter(
    !is.na(subject),
    is.finite(trial),
    is.finite(block),
    is.finite(arousal),
    is.finite(valence),
    is.finite(time_bin_s),
    is.finite(pupil_pct_change_sound_pre500to0_bin)
  )

# ============================================================
# 5. Check time-bin counts and create trial-level z-scored predictors
# ============================================================
trial_bin_counts <- M %>%
  group_by(subject, trial, block) %>%
  summarise(
    n_time_bins = n_distinct(time_bin_s),
    min_time_bin_s = min(time_bin_s, na.rm = TRUE),
    max_time_bin_s = max(time_bin_s, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  )

trial_bin_count_summary <- trial_bin_counts %>%
  summarise(
    n_trials = n(),
    min_n_bins = min(n_time_bins, na.rm = TRUE),
    max_n_bins = max(n_time_bins, na.rm = TRUE),
    mean_n_bins = mean(n_time_bins, na.rm = TRUE),
    sd_n_bins = sd(n_time_bins, na.rm = TRUE),
    n_unique_bin_counts = n_distinct(n_time_bins),
    min_start_time = min(min_time_bin_s, na.rm = TRUE),
    max_start_time = max(min_time_bin_s, na.rm = TRUE),
    min_end_time = min(max_time_bin_s, na.rm = TRUE),
    max_end_time = max(max_time_bin_s, na.rm = TRUE)
  )

trial_bin_count_distribution <- trial_bin_counts %>%
  count(n_time_bins, name = "n_trials") %>%
  arrange(n_time_bins)

mode_n_bins <- trial_bin_count_distribution %>%
  arrange(desc(n_trials), n_time_bins) %>%
  slice(1) %>%
  pull(n_time_bins)

trial_bin_count_problem_trials <- trial_bin_counts %>%
  filter(n_time_bins != mode_n_bins) %>%
  arrange(subject, trial)

write_csv(trial_bin_counts, file.path(OUTPUT_DIR, "trial_time_bin_counts.csv"))
write_csv(trial_bin_count_summary, file.path(OUTPUT_DIR, "trial_time_bin_count_summary.csv"))
write_csv(trial_bin_count_distribution, file.path(OUTPUT_DIR, "trial_time_bin_count_distribution.csv"))
write_csv(trial_bin_count_problem_trials, file.path(OUTPUT_DIR, "trial_time_bin_count_problem_trials.csv"))

cat("\nTime-bin count summary:\n")
print(trial_bin_count_summary)

if (trial_bin_count_summary$n_unique_bin_counts[1] == 1) {
  cat("All subject x trial units have the same number of time bins.\n")
} else {
  cat("WARNING: Some subject x trial units have different numbers of time bins.\n")
  cat("See trial_time_bin_count_problem_trials.csv\n")
}

trial_affect_z <- M %>%
  distinct(subject, trial, block, arousal, valence) %>%
  arrange(subject, trial) %>%
  group_by(subject) %>%
  mutate(
    arousal_z = safe_zscore(arousal),
    valence_z = safe_zscore(valence),
    arousal_subject_mean_raw = mean(arousal, na.rm = TRUE),
    arousal_subject_sd_raw = sd(arousal, na.rm = TRUE),
    valence_subject_mean_raw = mean(valence, na.rm = TRUE),
    valence_subject_sd_raw = sd(valence, na.rm = TRUE),
    n_trials_for_subject_z = sum(is.finite(arousal) & is.finite(valence))
  ) %>%
  ungroup()

trial_affect_z_summary <- trial_affect_z %>%
  summarise(
    zscoring = "within_participant_trial_level",
    n_trials = n(),
    n_subjects = n_distinct(subject),
    arousal_mean_raw = mean(arousal, na.rm = TRUE),
    arousal_sd_raw = sd(arousal, na.rm = TRUE),
    valence_mean_raw = mean(valence, na.rm = TRUE),
    valence_sd_raw = sd(valence, na.rm = TRUE),
    arousal_z_mean_grand = mean(arousal_z, na.rm = TRUE),
    arousal_z_sd_grand = sd(arousal_z, na.rm = TRUE),
    valence_z_mean_grand = mean(valence_z, na.rm = TRUE),
    valence_z_sd_grand = sd(valence_z, na.rm = TRUE),
    arousal_valence_r = cor(arousal_z, valence_z, use = "complete.obs")
  )

trial_affect_z_by_subject_summary <- trial_affect_z %>%
  group_by(subject) %>%
  summarise(
    n_trials = n(),
    arousal_mean_raw = mean(arousal, na.rm = TRUE),
    arousal_sd_raw = sd(arousal, na.rm = TRUE),
    valence_mean_raw = mean(valence, na.rm = TRUE),
    valence_sd_raw = sd(valence, na.rm = TRUE),
    arousal_z_mean = mean(arousal_z, na.rm = TRUE),
    arousal_z_sd = sd(arousal_z, na.rm = TRUE),
    valence_z_mean = mean(valence_z, na.rm = TRUE),
    valence_z_sd = sd(valence_z, na.rm = TRUE),
    arousal_valence_r = cor(arousal_z, valence_z, use = "complete.obs"),
    .groups = "drop"
  )

write_csv(trial_affect_z, file.path(OUTPUT_DIR, "trial_level_affect_zscores.csv"))
write_csv(trial_affect_z_summary, file.path(OUTPUT_DIR, "trial_level_affect_zscore_summary.csv"))
write_csv(trial_affect_z_by_subject_summary, file.path(OUTPUT_DIR, "trial_level_affect_zscore_by_subject_summary.csv"))

M <- M %>%
  select(-any_of(c("arousal_z", "valence_z"))) %>%
  left_join(
    trial_affect_z %>% select(subject, trial, arousal_z, valence_z),
    by = c("subject", "trial")
  ) %>%
  filter(is.finite(arousal_z), is.finite(valence_z))

# ============================================================
# 6. Grand-average pupil time course from the master table
# ============================================================
subj_pupil_main <- M %>%
  group_by(subject, block, time_bin_s) %>%
  summarise(
    pupil_mean = mean(pupil_pct_change_sound_pre500to0_bin, na.rm = TRUE),
    .groups = "drop"
  )

grand_pupil_main <- subj_pupil_main %>%
  group_by(block, time_bin_s) %>%
  summarise(
    mean = mean(pupil_mean, na.rm = TRUE),
    ci95 = ci95_t(pupil_mean),
    n_subjects = sum(is.finite(pupil_mean)),
    .groups = "drop"
  )

grand_pupil_main_wide <- grand_pupil_main %>%
  select(block, time_bin_s, mean, ci95, n_subjects) %>%
  mutate(block = paste0("block", block)) %>%
  pivot_wider(
    names_from = block,
    values_from = c(mean, ci95, n_subjects),
    names_sep = "_"
  ) %>%
  rename(trial_time = time_bin_s)

p_main_pupil <- plot_two_ci_lines(
  df = grand_pupil_main_wide,
  xcol = "trial_time",
  y1 = "mean_block1",
  ci1 = "ci95_block1",
  y2 = "mean_block2",
  ci2 = "ci95_block2",
  title_txt = "Grand-average pupil responses following sound onset",
  ylab_txt = "Pupil (% change from baseline)"
)

ggsave(
  filename = file.path(OUTPUT_DIR, "figure_main_pupil_timecourse.pdf"),
  plot = p_main_pupil,
  width = 8.5,
  height = 4.8
)

write_csv(grand_pupil_main_wide, file.path(OUTPUT_DIR, "grand_mean_pupil_main.csv"))

# ============================================================
# 7. Bin-wise mixed-effects analysis
# ============================================================
Tmix <- M %>%
  mutate(
    block_ec = case_when(
      block == 1 ~ -0.5,
      block == 2 ~ 0.5,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(is.finite(arousal_z), is.finite(valence_z), is.finite(block_ec))

Tmix$subject <- factor(Tmix$subject)
time_bins <- sort(unique(Tmix$time_bin_s))

result_list <- vector("list", length(time_bins))

for (i in seq_along(time_bins)) {
  tb <- time_bins[i]
  cat(sprintf("Processing bin %d / %d : t = %.3f s\n", i, length(time_bins), tb))

  D <- Tmix %>% filter(time_bin_s == tb)

  n_rows <- nrow(D)
  n_subjects <- n_distinct(D$subject)
  n_trials <- D %>% distinct(subject, trial) %>% nrow()

  one_row <- tibble(
    time_bin_s = as.numeric(tb),
    n_rows = as.integer(n_rows),
    n_subjects = as.integer(n_subjects),
    n_trials = as.integer(n_trials),
    model_converged = FALSE,
    note = "",

    beta_intercept = NA_real_,
    beta_arousal = NA_real_,
    beta_valence = NA_real_,
    beta_block = NA_real_,
    beta_arousal_block = NA_real_,
    beta_valence_block = NA_real_,

    se_intercept = NA_real_,
    se_arousal = NA_real_,
    se_valence = NA_real_,
    se_block = NA_real_,
    se_arousal_block = NA_real_,
    se_valence_block = NA_real_,

    df_intercept = NA_real_,
    df_arousal = NA_real_,
    df_valence = NA_real_,
    df_block = NA_real_,
    df_arousal_block = NA_real_,
    df_valence_block = NA_real_,

    t_intercept = NA_real_,
    t_arousal = NA_real_,
    t_valence = NA_real_,
    t_block = NA_real_,
    t_arousal_block = NA_real_,
    t_valence_block = NA_real_,

    p_intercept = NA_real_,
    p_arousal = NA_real_,
    p_valence = NA_real_,
    p_block = NA_real_,
    p_arousal_block = NA_real_,
    p_valence_block = NA_real_
  )

  if (n_rows < MIN_ROWS_PER_BIN || n_subjects < MIN_SUBJECTS_PER_BIN) {
    one_row$note <- "too_few_rows_or_subjects"
    result_list[[i]] <- one_row
    next
  }

  fit_warnings <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      lmer(formula = MODEL_FORMULA, data = D, REML = TRUE),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    one_row$model_converged <- FALSE
    one_row$note <- paste0("fit_error: ", conditionMessage(fit))
    result_list[[i]] <- one_row
    next
  }

  if (length(fit_warnings) > 0) {
    one_row$note <- paste0("fit_warning: ", paste(unique(fit_warnings), collapse = " | "))
  }

  coef_df <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) NULL)

  if (is.null(coef_df) || nrow(coef_df) == 0) {
    one_row$model_converged <- FALSE
    one_row$note <- "summary_coefficients_unavailable"
    result_list[[i]] <- one_row
    next
  }

  s_intercept <- get_coef_stats(coef_df, "(Intercept)")
  s_arousal <- get_coef_stats(coef_df, "arousal_z")
  s_valence <- get_coef_stats(coef_df, "valence_z")
  s_block <- get_coef_stats(coef_df, "block_ec")
  s_arousal_block <- get_coef_stats(coef_df, "arousal_z:block_ec")
  s_valence_block <- get_coef_stats(coef_df, "valence_z:block_ec")

  one_row$model_converged <- TRUE

  one_row$beta_intercept <- s_intercept$estimate
  one_row$beta_arousal <- s_arousal$estimate
  one_row$beta_valence <- s_valence$estimate
  one_row$beta_block <- s_block$estimate
  one_row$beta_arousal_block <- s_arousal_block$estimate
  one_row$beta_valence_block <- s_valence_block$estimate

  one_row$se_intercept <- s_intercept$se
  one_row$se_arousal <- s_arousal$se
  one_row$se_valence <- s_valence$se
  one_row$se_block <- s_block$se
  one_row$se_arousal_block <- s_arousal_block$se
  one_row$se_valence_block <- s_valence_block$se

  one_row$df_intercept <- s_intercept$df
  one_row$df_arousal <- s_arousal$df
  one_row$df_valence <- s_valence$df
  one_row$df_block <- s_block$df
  one_row$df_arousal_block <- s_arousal_block$df
  one_row$df_valence_block <- s_valence_block$df

  one_row$t_intercept <- s_intercept$t
  one_row$t_arousal <- s_arousal$t
  one_row$t_valence <- s_valence$t
  one_row$t_block <- s_block$t
  one_row$t_arousal_block <- s_arousal_block$t
  one_row$t_valence_block <- s_valence_block$t

  one_row$p_intercept <- s_intercept$p
  one_row$p_arousal <- s_arousal$p
  one_row$p_valence <- s_valence$p
  one_row$p_block <- s_block$p
  one_row$p_arousal_block <- s_arousal_block$p
  one_row$p_valence_block <- s_valence_block$p

  result_list[[i]] <- one_row
}

Rmix <- bind_rows(result_list) %>%
  mutate(
    arousal_slope_block1 = beta_arousal - 0.5 * beta_arousal_block,
    arousal_slope_block2 = beta_arousal + 0.5 * beta_arousal_block,
    valence_slope_block1 = beta_valence - 0.5 * beta_valence_block,
    valence_slope_block2 = beta_valence + 0.5 * beta_valence_block,
    block2_minus_block1_arousal_slope = beta_arousal_block,
    block2_minus_block1_valence_slope = beta_valence_block,

    direction_arousal_main = vapply(beta_arousal, interpret_direction, character(1)),
    direction_valence_main = vapply(beta_valence, interpret_direction, character(1)),
    direction_block_main = vapply(beta_block, interpret_direction, character(1)),
    direction_arousal_block_diff = vapply(beta_arousal_block, interpret_direction, character(1)),
    direction_valence_block_diff = vapply(beta_valence_block, interpret_direction, character(1)),

    q_intercept = p.adjust(p_intercept, method = "BH"),
    q_arousal = p.adjust(p_arousal, method = "BH"),
    q_valence = p.adjust(p_valence, method = "BH"),
    q_block = p.adjust(p_block, method = "BH"),
    q_arousal_block = p.adjust(p_arousal_block, method = "BH"),
    q_valence_block = p.adjust(p_valence_block, method = "BH"),

    sig_intercept_q05 = q_intercept < ALPHA_Q,
    sig_arousal_q05 = q_arousal < ALPHA_Q,
    sig_valence_q05 = q_valence < ALPHA_Q,
    sig_block_q05 = q_block < ALPHA_Q,
    sig_arousal_block_q05 = q_arousal_block < ALPHA_Q,
    sig_valence_block_q05 = q_valence_block < ALPHA_Q
  )

Rmix$effect_pattern_q05 <- mapply(
  summarize_effect_pattern,
  Rmix$sig_arousal_q05,
  Rmix$sig_valence_q05,
  Rmix$sig_block_q05,
  Rmix$sig_arousal_block_q05,
  Rmix$sig_valence_block_q05
)

write_csv(Rmix, file.path(OUTPUT_DIR, "binwise_mixed_model_results_raw.csv"))
write_csv(Rmix, file.path(OUTPUT_DIR, "binwise_mixed_model_results_fdr.csv"))

# ============================================================
# 8. Mixed-model beta figures
# ============================================================
p_beta_1 <- plot_beta_panel(Rmix, "beta_arousal", "q_arousal", "Arousal", "beta(arousal)")
p_beta_2 <- plot_beta_panel(Rmix, "beta_valence", "q_valence", "Valence", "beta(valence)")
p_beta_3 <- plot_beta_panel(Rmix, "beta_block", "q_block", "Block", "beta(block)")
p_beta_4 <- plot_beta_panel(Rmix, "beta_arousal_block", "q_arousal_block", "Arousal x block", "beta(arousal x block)")
p_beta_5 <- plot_beta_panel(Rmix, "beta_valence_block", "q_valence_block", "Valence x block", "beta(valence x block)")

ggsave(file.path(OUTPUT_DIR, "figure_main_mixed_beta_timecourses_1_arousal.pdf"), p_beta_1, width = 8.5, height = 3.6)
ggsave(file.path(OUTPUT_DIR, "figure_main_mixed_beta_timecourses_2_valence.pdf"), p_beta_2, width = 8.5, height = 3.6)
ggsave(file.path(OUTPUT_DIR, "figure_main_mixed_beta_timecourses_3_block.pdf"), p_beta_3, width = 8.5, height = 3.6)
ggsave(file.path(OUTPUT_DIR, "figure_main_mixed_beta_timecourses_4_arousal_block.pdf"), p_beta_4, width = 8.5, height = 3.6)
ggsave(file.path(OUTPUT_DIR, "figure_main_mixed_beta_timecourses_5_valence_block.pdf"), p_beta_5, width = 8.5, height = 3.6)

# Four-panel manuscript-style figure: arousal, valence, arousal x block, valence x block.
save_four_panel_pdf(
  p_beta_1,
  p_beta_2,
  p_beta_4,
  p_beta_5,
  file.path(OUTPUT_DIR, "figure_main_mixed_beta_timecourses_4panel.pdf")
)

# ============================================================
# 9. Save run info
# ============================================================
run_info_file <- file.path(OUTPUT_DIR, "run_info.txt")
con <- file(run_info_file, open = "wt", encoding = "UTF-8")

writeLines(sprintf("MASTER_FILE=%s", MASTER_FILE), con)
writeLines(sprintf("OUTPUT_DIR=%s", OUTPUT_DIR), con)
writeLines(sprintf("DEPENDENT_VAR=%s", DEPENDENT_VAR), con)
writeLines(sprintf("MIN_ROWS_PER_BIN=%d", MIN_ROWS_PER_BIN), con)
writeLines(sprintf("MIN_SUBJECTS_PER_BIN=%d", MIN_SUBJECTS_PER_BIN), con)
writeLines("MIXED_MODEL_ENGINE=lmerTest::lmer", con)
writeLines(sprintf("MIXED_MODEL_FORMULA=%s", deparse(MODEL_FORMULA)), con)
writeLines("BLOCK_CODING=effect_coding", con)
writeLines("BLOCK1_CODE=-0.5", con)
writeLines("BLOCK2_CODE=+0.5", con)
writeLines("P_VALUES=t_distribution_based_from_lmerTest", con)
writeLines("FDR=Benjamini-Hochberg separately by effect across time bins", con)
writeLines("AROUSAL_VALENCE_ZSCORING=within-participant trial-level z-score", con)
writeLines("AROUSAL_Z_SOURCE=computed within each subject after reducing to distinct subject x trial rows", con)
writeLines("VALENCE_Z_SOURCE=computed within each subject after reducing to distinct subject x trial rows", con)

writeLines(sprintf("TIME_BIN_COUNT_N_TRIALS=%d", trial_bin_count_summary$n_trials[1]), con)
writeLines(sprintf("TIME_BIN_COUNT_MIN=%d", trial_bin_count_summary$min_n_bins[1]), con)
writeLines(sprintf("TIME_BIN_COUNT_MAX=%d", trial_bin_count_summary$max_n_bins[1]), con)
writeLines(sprintf("TIME_BIN_COUNT_UNIQUE=%d", trial_bin_count_summary$n_unique_bin_counts[1]), con)
writeLines(sprintf(
  "TIME_BIN_START_RANGE=%.6f_to_%.6f",
  trial_bin_count_summary$min_start_time[1],
  trial_bin_count_summary$max_start_time[1]
), con)
writeLines(sprintf(
  "TIME_BIN_END_RANGE=%.6f_to_%.6f",
  trial_bin_count_summary$min_end_time[1],
  trial_bin_count_summary$max_end_time[1]
), con)
writeLines(sprintf(
  "AROUSAL_VALENCE_CORRELATION_TRIAL_LEVEL=%.6f",
  trial_affect_z_summary$arousal_valence_r[1]
), con)
writeLines("WAVEFORM_CI=95% t-based CI across subjects", con)

close(con)

cat("\nSaved figures and tables to:\n", OUTPUT_DIR, "\n")
cat("=== Main bin-wise analysis complete ===\n")
