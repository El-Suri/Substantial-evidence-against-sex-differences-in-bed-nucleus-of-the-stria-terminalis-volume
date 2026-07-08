args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/02_create_prior_table.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("dplyr", "metafor"))

prior_dir <- file.path(root, "data", "prior_sources")
results_dir <- file.path(root, "results")
figures_dir <- file.path(root, "figures")
ensure_dir(results_dir)
ensure_dir(figures_dir)

n_boot <- as.integer(read_arg("n_boot", "5000"))
guma_participant_file <- read_arg("guma_participant_file")
guma_adjusted_file <- read_arg("guma_adjusted_effect", file.path(prior_dir, "Guma_2024_adjusted_effect.csv"))
set.seed(123)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------
# The prior is built from standardized male-female BNST differences. Where
# participant-level data are available, we estimate covariate-adjusted effects
# using a linear model and standardize the adjusted contrast by the model's
# residual SD. Where only paper summaries are available, we derive Cohen's d
# from the reported means, SDs/SEMs, and sample sizes.

calc_raw_d <- function(data, outcome, group, male_level, female_level) {
  y_m <- data[[outcome]][data[[group]] == male_level]
  y_f <- data[[outcome]][data[[group]] == female_level]
  n_m <- sum(!is.na(y_m))
  n_f <- sum(!is.na(y_f))
  sd_m <- stats::sd(y_m, na.rm = TRUE)
  sd_f <- stats::sd(y_f, na.rm = TRUE)
  sd_pooled <- sqrt(((n_m - 1) * sd_m^2 + (n_f - 1) * sd_f^2) / (n_m + n_f - 2))
  d <- (mean(y_m, na.rm = TRUE) - mean(y_f, na.rm = TRUE)) / sd_pooled
  se_d <- sqrt((n_m + n_f) / (n_m * n_f) + d^2 / (2 * (n_m + n_f - 2)))

  data.frame(
    n_male = n_m,
    n_female = n_f,
    effect_size = d,
    effect_se = se_d
  )
}

calc_smd_from_summary <- function(mean_m, mean_f, sd_m, sd_f, n_m, n_f) {
  sd_pooled <- sqrt(((n_m - 1) * sd_m^2 + (n_f - 1) * sd_f^2) / (n_m + n_f - 2))
  d <- (mean_m - mean_f) / sd_pooled
  se_d <- sqrt((n_m + n_f) / (n_m * n_f) + d^2 / (2 * (n_m + n_f - 2)))
  data.frame(
    n_male = n_m,
    n_female = n_f,
    effect_size = d,
    effect_se = se_d
  )
}

combine_subgroups <- function(n1, mean1, sd1, n2, mean2, sd2) {
  pooled_mean <- (n1 * mean1 + n2 * mean2) / (n1 + n2)
  pooled_var <- (
    ((n1 - 1) * sd1^2) +
      ((n2 - 1) * sd2^2) +
      n1 * (mean1 - pooled_mean)^2 +
      n2 * (mean2 - pooled_mean)^2
  ) / (n1 + n2 - 1)
  data.frame(n = n1 + n2, mean = pooled_mean, sd = sqrt(pooled_var))
}

adjusted_d_from_lm <- function(data, outcome, group, covariates,
                               male_level, female_level, seed, study_label) {
  vars <- c(outcome, group, covariates)
  df <- data[stats::complete.cases(data[, vars]), vars, drop = FALSE]
  df[[group]] <- stats::relevel(factor(df[[group]]), ref = female_level)

  fmla <- stats::as.formula(paste(outcome, "~", paste(c(group, covariates), collapse = " + ")))
  fit <- stats::lm(fmla, data = df)

  # Predict male and female values at the sample mean of each covariate.
  new_f <- data.frame(row.names = 1)
  new_m <- data.frame(row.names = 1)
  new_f[[group]] <- factor(female_level, levels = levels(df[[group]]))
  new_m[[group]] <- factor(male_level, levels = levels(df[[group]]))
  for (v in covariates) {
    new_f[[v]] <- mean(df[[v]], na.rm = TRUE)
    new_m[[v]] <- mean(df[[v]], na.rm = TRUE)
  }

  adjusted_mean_male <- as.numeric(stats::predict(fit, new_m))
  adjusted_mean_female <- as.numeric(stats::predict(fit, new_f))
  diff_adj <- adjusted_mean_male - adjusted_mean_female
  d_adj <- diff_adj / stats::sigma(fit)

  # Nonparametric bootstrap SE for the adjusted standardized contrast.
  set.seed(seed)
  boot_vals <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
    df_b <- df[idx, , drop = FALSE]
    df_b[[group]] <- factor(df_b[[group]], levels = levels(df[[group]]))
    if (length(unique(df_b[[group]])) < 2) next
    fit_b <- try(stats::lm(fmla, data = df_b), silent = TRUE)
    if (inherits(fit_b, "try-error")) next
    boot_vals[b] <- as.numeric(stats::predict(fit_b, new_m) - stats::predict(fit_b, new_f)) / stats::sigma(fit_b)
  }
  boot_vals <- boot_vals[is.finite(boot_vals)]

  data.frame(
    study = study_label,
    n_total = nrow(df),
    n_male = sum(df[[group]] == male_level),
    n_female = sum(df[[group]] == female_level),
    formula = deparse(fmla),
    adjusted_mean_male = adjusted_mean_male,
    adjusted_mean_female = adjusted_mean_female,
    adjusted_difference = diff_adj,
    residual_sd = stats::sigma(fit),
    effect_size = d_adj,
    effect_se = stats::sd(boot_vals),
    n_boot_valid = length(boot_vals),
    ci_lower = stats::quantile(boot_vals, 0.025, na.rm = TRUE),
    ci_upper = stats::quantile(boot_vals, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

add_result <- function(study, source_type, sample_restriction, effect_method,
                       covariates, n_male, n_female, effect_size, effect_se, notes) {
  data.frame(
    study = study,
    source_type = source_type,
    sample_restriction = sample_restriction,
    effect_method = effect_method,
    covariates = covariates,
    n_male = n_male,
    n_female = n_female,
    effect_size = effect_size,
    effect_se = effect_se,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

require_file <- function(path, description) {
  if (!file.exists(path)) {
    stop(
      "Missing ", description, ": ", path, "\n",
      "See data/prior_sources/README.md for required prior-source inputs.",
      call. = FALSE
    )
  }
}

# -------------------------------------------------------------------------
# Participant-level re-analyses
# -------------------------------------------------------------------------

allen_path <- file.path(prior_dir, "Allen_Gorski.csv")
require_file(allen_path, "Allen & Gorski participant-level table")
allen <- read.csv(allen_path, check.names = TRUE)
allen <- allen[allen$Age >= 18, , drop = FALSE]
allen_res <- adjusted_d_from_lm(
  data = allen,
  outcome = "BNST.Volume",
  group = "Sex",
  covariates = c("Brain.Weight", "Age"),
  male_level = "m",
  female_level = "f",
  seed = 101,
  study_label = "Allen & Gorski (1990)"
)

chung_path <- file.path(prior_dir, "Chung_etal.csv")
require_file(chung_path, "Chung participant-level table")
chung <- read.csv(chung_path, check.names = TRUE)
chung$Age <- as.numeric(gsub(" years", "", chung$Age))
chung <- chung[chung$Age >= 18, , drop = FALSE]
chung_res <- adjusted_d_from_lm(
  data = chung,
  outcome = "BNST.Volume",
  group = "Sex",
  covariates = c("Brain.Weight", "Age"),
  male_level = "m",
  female_level = "f",
  seed = 102,
  study_label = "Chung et al. (2002)"
)

# Guma et al. participant-level HCP data cannot be included publicly. If a
# local RDS is supplied, recompute the adjusted effect. Otherwise use a cached
# study-level adjusted-effect CSV that contains only n/effect/SE.
if (!is.null(guma_participant_file) && file.exists(guma_participant_file)) {
  guma <- readRDS(guma_participant_file)
  guma <- subset(guma, euler > -200 & QC_include == 1)
  guma$Sex <- factor(guma$Sex, levels = c("F", "M"))
  guma$BNST_bilateral <- guma$L_BNST + guma$R_BNST
  guma_res <- adjusted_d_from_lm(
    data = guma,
    outcome = "BNST_bilateral",
    group = "Sex",
    covariates = c("AGE_cent", "BrainSegVolNotVent.y", "euler"),
    male_level = "M",
    female_level = "F",
    seed = 103,
    study_label = "Guma et al. (2024)"
  )
} else {
  require_file(guma_adjusted_file, "Guma adjusted-effect table")
  guma_res <- read.csv(guma_adjusted_file, check.names = TRUE)
}

# -------------------------------------------------------------------------
# Published summary-statistic derivations
# -------------------------------------------------------------------------

# Zhou et al. report mean +/- SEM. We combine cisgender heterosexual and
# homosexual male groups, exclude hormone-treated transgender participants, and
# compare the combined male group with female participants.
zhou_male <- combine_subgroups(
  n1 = 15, mean1 = 2.49, sd1 = 0.16 * sqrt(15),
  n2 = 9, mean2 = 2.81, sd2 = 0.20 * sqrt(9)
)
zhou <- calc_smd_from_summary(
  mean_m = zhou_male$mean,
  mean_f = 1.73,
  sd_m = zhou_male$sd,
  sd_f = 0.13 * sqrt(12),
  n_m = zhou_male$n,
  n_f = 12
)

slabe <- calc_smd_from_summary(
  mean_m = 4.72,
  mean_f = 4.27,
  sd_m = 1.56,
  sd_f = 1.57,
  n_m = 22,
  n_f = 12
)

# Neudorfer et al. did not report sex-specific sample sizes. The manuscript
# estimates the split from the broader HCP Young Adult sex proportions
# described in Guma et al. (597 female, 496 male).
neud_total_n <- 990
neud_n_f <- round((597 / (597 + 496)) * neud_total_n)
neud_n_m <- neud_total_n - neud_n_f
neud <- calc_smd_from_summary(
  mean_m = 8.02e-05,
  mean_f = 7.97e-05,
  sd_m = 9.13e-06,
  sd_f = 9.32e-06,
  n_m = neud_n_m,
  n_f = neud_n_f
)

# -------------------------------------------------------------------------
# Combine effect-size table and fit the random-effects prior meta-analysis.
# -------------------------------------------------------------------------

prior_table <- rbind(
  add_result("Allen & Gorski (1990)", "participant-level re-analysis", "Adults only (>=18 years)", "Adjusted standardized mean difference from lm", "Brain.Weight, Age", allen_res$n_male, allen_res$n_female, allen_res$effect_size, allen_res$effect_se, "Effect derived from regression-adjusted male-female contrast divided by residual SD."),
  add_result("Chung et al. (2002)", "participant-level re-analysis", "Adults only (>=18 years)", "Adjusted standardized mean difference from lm", "Brain.Weight, Age", chung_res$n_male, chung_res$n_female, chung_res$effect_size, chung_res$effect_se, "Effect derived from regression-adjusted male-female contrast divided by residual SD."),
  add_result("Guma et al. (2024)", "participant-level re-analysis", "Paper-matched HCP QC subset", "Adjusted standardized mean difference from lm", "AGE_cent, BrainSegVolNotVent.y, euler", guma_res$n_male[1], guma_res$n_female[1], guma_res$effect_size[1], guma_res$effect_se[1], "Effect derived from regression-adjusted male-female contrast divided by residual SD."),
  add_result("Zhou et al. (1995)", "published summary statistics", "Adult cisgender participants; transgender hormone-treated participants excluded", "Cohen's d from reported means, SEMs, and reconstructed subgroup SDs", "None available", zhou$n_male, zhou$n_female, zhou$effect_size, zhou$effect_se, "Male group combines heterosexual and homosexual cisgender men; SEMs converted to SDs before pooling."),
  add_result("Slabe et al. (2023)", "published summary statistics", "As reported in paper", "Cohen's d from reported means, SDs, and n", "No global size adjustment available", slabe$n_male, slabe$n_female, slabe$effect_size, slabe$effect_se, "No TBV/brain-weight-adjusted effect available."),
  add_result("Neudorfer et al. (2020)", "published TBV-normalised summary statistics", "As reported in paper", "Cohen's d from TBV-normalised means, SDs, and estimated sex-specific n", "TBV-normalised values reported by paper", neud$n_male, neud$n_female, neud$effect_size, neud$effect_se, "Sex split estimated from Guma et al. HCP sex proportions because Neudorfer did not report sex-specific sample counts.")
)

write.csv(prior_table, file.path(results_dir, "BNST_prior_effect_sizes_all_studies.csv"), row.names = FALSE)
write.csv(allen_res, file.path(results_dir, "Allen_Gorski_1990_adjusted_effect.csv"), row.names = FALSE)
write.csv(chung_res, file.path(results_dir, "Chung_2002_adjusted_effect.csv"), row.names = FALSE)
write.csv(guma_res, file.path(results_dir, "Guma_2024_adjusted_effect.csv"), row.names = FALSE)

meta_input <- prior_table[, c("study", "effect_size", "effect_se")]
names(meta_input) <- c("study", "yi", "sei")
meta_input <- meta_input[is.finite(meta_input$yi) & is.finite(meta_input$sei), ]
write.csv(meta_input, file.path(results_dir, "BNST_prior_meta_input.csv"), row.names = FALSE)

primary <- metafor::rma(yi = yi, sei = sei, data = meta_input, method = "REML")

summary_out <- data.frame(
  analysis = "primary",
  pooled_effect = as.numeric(primary$b),
  pooled_se = primary$se,
  ci_lower = primary$ci.lb,
  ci_upper = primary$ci.ub,
  tau2 = primary$tau2,
  I2 = primary$I2,
  Q = primary$QE,
  Q_p = primary$QEp
)

meta_no_neud <- meta_input[!grepl("Neudorfer", meta_input$study, ignore.case = TRUE), ]
meta_no_zhou <- meta_input[!grepl("Zhou", meta_input$study, ignore.case = TRUE), ]
meta_adjusted_only <- meta_input[meta_input$study %in% c("Allen & Gorski (1990)", "Chung et al. (2002)", "Guma et al. (2024)"), ]

sensitivity <- list(
  exclude_neudorfer = metafor::rma(yi = yi, sei = sei, data = meta_no_neud, method = "REML"),
  exclude_zhou = metafor::rma(yi = yi, sei = sei, data = meta_no_zhou, method = "REML"),
  adjusted_only = metafor::rma(yi = yi, sei = sei, data = meta_adjusted_only, method = "REML")
)

sensitivity_table <- do.call(rbind, lapply(names(sensitivity), function(label) {
  fit <- sensitivity[[label]]
  data.frame(
    analysis = label,
    pooled_effect = as.numeric(fit$b),
    pooled_se = fit$se,
    ci_lower = fit$ci.lb,
    ci_upper = fit$ci.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    Q = fit$QE,
    Q_p = fit$QEp
  )
}))
write.csv(rbind(summary_out, sensitivity_table), file.path(results_dir, "BNST_prior_meta_summary.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Prior/meta-analysis figures used for checking and manuscript drafting.
# -------------------------------------------------------------------------

grDevices::pdf(file.path(figures_dir, "BNST_meta_analysis_forest_plot.pdf"), width = 6, height = 5)
metafor::forest(primary, slab = meta_input$study, xlab = "Standardized mean difference", mlab = "Random-effects model")
grDevices::dev.off()

grDevices::pdf(file.path(figures_dir, "BNST_meta_analysis_forest_plots_sensitivity.pdf"), width = 6, height = 5)
metafor::forest(sensitivity$exclude_neudorfer, slab = meta_no_neud$study, xlab = "Standardized mean difference", mlab = "Exclude Neudorfer")
metafor::forest(sensitivity$exclude_zhou, slab = meta_no_zhou$study, xlab = "Standardized mean difference", mlab = "Exclude Zhou")
metafor::forest(sensitivity$adjusted_only, slab = meta_adjusted_only$study, xlab = "Standardized mean difference", mlab = "Adjusted-only studies")
grDevices::dev.off()

# influence() is an S3 generic; the rma.uni method is registered by metafor
# but is not exported as metafor::influence().
inf <- stats::influence(primary)
grDevices::pdf(file.path(figures_dir, "BNST_meta_analysis_influence_plot.pdf"), width = 7, height = 5)
plot(inf)
grDevices::dev.off()

capture.output(
  {
    cat("Primary random-effects meta-analysis\n")
    print(summary(primary))
    cat("\nSensitivity analyses\n")
    print(rbind(summary_out, sensitivity_table))
    cat("\nInfluence diagnostics\n")
    print(inf)
  },
  file = file.path(results_dir, "BNST_meta_analysis_summary.txt")
)

message("Wrote prior effect table, meta-analysis summaries, and prior figures.")
message("Primary pooled d = ", round(as.numeric(primary$b), 4), "; SE = ", round(primary$se, 4))
