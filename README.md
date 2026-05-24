# Main bin-wise pupil mixed-effects analysis

This repository contains the filtered master table and R script used to reproduce the main time-resolved pupil analysis for the manuscript:

**Emotional Arousal Predicts Pupil Responses When Attended Luminance Changes**

The analysis tests whether trial-by-trial subjective arousal and valence ratings predict post-sound pupil responses while participants shift attention between bright and dark moving dot fields.

## Repository contents

Recommended layout:

```text
repository/
  README.md
  data/
    all_subjects_master_table_main_analysis_filtered.csv
  scripts/
    05_main_binwise_mixed_model_from_master_csv.R
  analysis_outputs/
```

The script can also be run from other layouts by passing the input CSV and output directory as command-line arguments.

## Input data

The script uses one input file only:

```text
all_subjects_master_table_main_analysis_filtered.csv
```

This file is the filtered master table generated after pupil preprocessing and trial-level exclusion. It contains 100-ms binned post-sound pupil responses and trial-level emotional ratings.

Required columns are:

```text
subject
trial
block
arousal
valence
time_bin_s
pupil_pct_change_sound_pre500to0_bin
```

where `pupil_pct_change_sound_pre500to0_bin` is the binned pupil response expressed as percentage change from the 500-ms pre-sound baseline.

## R requirements

The script requires R and the following packages:

```r
install.packages(c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "lmerTest"
))
```

The `grid` package is also used and is included with base R.

## How to run

From the repository root, run:

```bash
Rscript scripts/05_main_binwise_mixed_model_from_master_csv.R \
  data/all_subjects_master_table_main_analysis_filtered.csv \
  analysis_outputs/main_binwise_mixed_model
```

If no command-line arguments are supplied, the script searches for the input file in the following locations:

```text
data/all_subjects_master_table_main_analysis_filtered.csv
all_subjects_master_table_main_analysis_filtered.csv
```

and writes results to:

```text
analysis_outputs/main_binwise_mixed_model/run_YYYYMMDD_HHMMSS
```

## Analysis summary

The script performs the main bin-wise mixed-effects analysis described in the manuscript.

1. Arousal and valence ratings are reduced to one value per subject and trial.
2. Arousal and valence are z-standardized within participant.
3. Block is effect-coded as:

```text
block 1 = -0.5
block 2 = +0.5
```

4. For each 100-ms time bin from 0 to 6 s after sound onset, the following mixed-effects model is fitted using `lmerTest::lmer`:

```text
pupil_pct_change_sound_pre500to0_bin ~
  arousal_z + valence_z + block_ec +
  arousal_z:block_ec + valence_z:block_ec +
  (1 | subject)
```

5. Time bins are modeled only when they include at least 10 observations and data from at least two subjects.
6. P values are adjusted across time bins using the Benjamini-Hochberg false discovery rate procedure separately for each effect.

## Outputs

Each run creates a timestamped output directory containing figures, result tables, and run metadata.

Main figures:

```text
figure_main_pupil_timecourse.pdf
figure_main_mixed_beta_timecourses_4panel.pdf
figure_main_mixed_beta_timecourses_1_arousal.pdf
figure_main_mixed_beta_timecourses_2_valence.pdf
figure_main_mixed_beta_timecourses_3_block.pdf
figure_main_mixed_beta_timecourses_4_arousal_block.pdf
figure_main_mixed_beta_timecourses_5_valence_block.pdf
```

Main result tables:

```text
grand_mean_pupil_main.csv
binwise_mixed_model_results_raw.csv
binwise_mixed_model_results_fdr.csv
trial_level_affect_zscores.csv
trial_level_affect_zscore_summary.csv
trial_level_affect_zscore_by_subject_summary.csv
trial_time_bin_counts.csv
trial_time_bin_count_summary.csv
trial_time_bin_count_distribution.csv
trial_time_bin_count_problem_trials.csv
run_info.txt
```

`run_info.txt` records the input file, output directory, model formula, block coding, FDR procedure, z-scoring procedure, and basic time-bin diagnostics.

## Scope of this script

This script is intentionally limited to analyses that can be reproduced from `all_subjects_master_table_main_analysis_filtered.csv` alone. It reproduces the main post-sound pupil time-course figure and the main time-resolved mixed-effects coefficient figures and tables.

It does not reproduce supplementary long-window pupil or OKN visualizations, because those require separate full time-series and OKN-derived files. It also does not reproduce the separate robustness-analysis scripts unless those scripts and their required input files are included separately in the repository.

## Notes for manuscript reporting

A concise description of this repository can be used in the manuscript as follows:

```text
The filtered master table used for the main bin-wise mixed-effects analysis and the corresponding R analysis script are available in the repository. The script uses all_subjects_master_table_main_analysis_filtered.csv as its only input, computes within-participant trial-level z-scored arousal and valence predictors, fits the time-resolved mixed-effects models, applies Benjamini-Hochberg FDR correction across time bins separately for each effect, and outputs the main pupil time-course figure, mixed-model coefficient time courses, and model-result tables.
```
