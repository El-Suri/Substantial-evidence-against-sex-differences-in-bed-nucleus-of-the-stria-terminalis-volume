args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/02_create_prior_table.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("metafor"))

n_boot <- as.integer(read_arg("n_boot", "5000"))
set.seed(123)

calc_smd_from_summary <- function(mean_m, mean_f, sd_m, sd_f, n_m, n_f) {
  sd_pooled <- sqrt(((n_m - 1) * sd_m^2 + (n_f - 1) * sd_f^2) / (n_m + n_f - 2))
  d <- (mean_m - mean_f) / sd_pooled
  se_d <- sqrt((n_m + n_f) / (n_m * n_f) + d^2 / (2 * (n_m + n_f - 2)))
  list(effect_size = d, effect_se = se_d)
}

combine_subgroups <- function(n1, mean1, sd1, n2, mean2, sd2) {
  pooled_mean <- (n1 * mean1 + n2 * mean2) / (n1 + n2)
  pooled_var <- (
    ((n1 - 1) * sd1^2) +
      ((n2 - 1) * sd2^2) +
      n1 * (mean1 - pooled_mean)^2 +
      n2 * (mean2 - pooled_mean)^2
  ) / (n1 + n2 - 1)
  list(n = n1 + n2, mean = pooled_mean, sd = sqrt(pooled_var))
}

adjusted_d_from_lm <- function(data, outcome, group, covariates, male_level, female_level, seed) {
  vars <- c(outcome, group, covariates)
  df <- data[complete.cases(data[, vars]), vars, drop = FALSE]
  df[[group]] <- stats::relevel(factor(df[[group]]), ref = female_level)
  fmla <- stats::as.formula(paste(outcome, "~", paste(c(group, covariates), collapse = " + ")))
  fit <- stats::lm(fmla, data = df)

  new_f <- data.frame(row.names = 1)
  new_m <- data.frame(row.names = 1)
  new_f[[group]] <- factor(female_level, levels = levels(df[[group]]))
  new_m[[group]] <- factor(male_level, levels = levels(df[[group]]))
  for (v in covariates) {
    new_f[[v]] <- mean(df[[v]], na.rm = TRUE)
    new_m[[v]] <- mean(df[[v]], na.rm = TRUE)
  }

  diff_adj <- as.numeric(stats::predict(fit, new_m) - stats::predict(fit, new_f))
  d_adj <- diff_adj / stats::sigma(fit)

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

  list(
    n_total = nrow(df),
    n_male = sum(df[[group]] == male_level),
    n_female = sum(df[[group]] == female_level),
    effect_size = d_adj,
    effect_se = stats::sd(boot_vals)
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

prior_dir <- file.path(root, "data", "prior_sources")

allen <- read.csv(file.path(prior_dir, "Allen_Gorski.csv"), check.names = TRUE)
allen <- allen[allen$Age >= 18, ]
allen_res <- adjusted_d_from_lm(allen, "BNST.Volume", "Sex", c("Brain.Weight", "Age"), "m", "f", 101)

chung <- read.csv(file.path(prior_dir, "Chung_etal.csv"), check.names = TRUE)
chung$Age <- as.numeric(gsub(" years", "", chung$Age))
chung <- chung[chung$Age >= 18, ]
chung_res <- adjusted_d_from_lm(chung, "BNST.Volume", "Sex", c("Brain.Weight", "Age"), "m", "f", 102)

guma <- read.csv(file.path(prior_dir, "Guma_2024_adjusted_effect.csv"), check.names = TRUE)

zhou_m <- combine_subgroups(
  n1 = 15, mean1 = 2.49, sd1 = 0.16 * sqrt(15),
  n2 = 9, mean2 = 2.81, sd2 = 0.20 * sqrt(9)
)
zhou <- calc_smd_from_summary(zhou_m$mean, 1.73, zhou_m$sd, 0.13 * sqrt(12), zhou_m$n, 12)

slabe <- calc_smd_from_summary(4.72, 4.27, 1.56, 1.57, 22, 12)

neud_n_total <- 990
neud_n_f <- round((597 / (597 + 496)) * neud_n_total)
neud_n_m <- neud_n_total - neud_n_f
neud <- calc_smd_from_summary(8.02e-05, 7.97e-05, 9.13e-06, 9.32e-06, neud_n_m, neud_n_f)

prior_table <- rbind(
  add_result("Allen & Gorski (1990)", "participant-level re-analysis", "Adults only (>=18 years)", "Adjusted standardized mean difference from lm", "Brain.Weight, Age", allen_res$n_male, allen_res$n_female, allen_res$effect_size, allen_res$effect_se, "Recomputed from included participant-level table."),
  add_result("Chung et al. (2002)", "participant-level re-analysis", "Adults only (>=18 years)", "Adjusted standardized mean difference from lm", "Brain.Weight, Age", chung_res$n_male, chung_res$n_female, chung_res$effect_size, chung_res$effect_se, "Recomputed from included participant-level table."),
  add_result("Guma et al. (2024)", "participant-level re-analysis", "Paper-matched HCP QC subset", "Adjusted standardized mean difference from lm", "AGE_cent, BrainSegVolNotVent.y, euler", guma$n_male[1], guma$n_female[1], guma$effect_size[1], guma$effect_se[1], "Cached adjusted effect; participant-level HCP RDS not redistributed."),
  add_result("Zhou et al. (1995)", "published summary statistics", "Adult cisgender participants; transgender hormone-treated participants excluded", "Cohen's d from reported means, SEMs, and reconstructed subgroup SDs", "None available", zhou_m$n, 12, zhou$effect_size, zhou$effect_se, "Male group combines heterosexual and homosexual cisgender men."),
  add_result("Slabe et al. (2023)", "published summary statistics", "As reported in paper", "Cohen's d from reported means, SDs, and n", "No global size adjustment available", 22, 12, slabe$effect_size, slabe$effect_se, "No TBV/brain-weight-adjusted effect available."),
  add_result("Neudorfer et al. (2020)", "published TBV-normalised summary statistics", "As reported in paper", "Cohen's d from TBV-normalised means, SDs, and estimated sex-specific n", "TBV-normalised values reported by paper", neud_n_m, neud_n_f, neud$effect_size, neud$effect_se, "Sex split estimated from HCP Young Adult open demographic proportions.")
)

results_dir <- file.path(root, "results")
ensure_dir(results_dir)
write.csv(prior_table, file.path(results_dir, "BNST_prior_effect_sizes_all_studies_recomputed.csv"), row.names = FALSE)

meta_input <- prior_table[, c("study", "effect_size", "effect_se")]
names(meta_input) <- c("study", "yi", "sei")
res <- metafor::rma(yi = yi, sei = sei, data = meta_input, method = "REML")

summary_out <- data.frame(
  pooled_effect = as.numeric(res$b),
  pooled_se = res$se,
  ci_lower = res$ci.lb,
  ci_upper = res$ci.ub,
  tau2 = res$tau2,
  I2 = res$I2
)
write.csv(summary_out, file.path(results_dir, "BNST_prior_meta_summary_recomputed.csv"), row.names = FALSE)

grDevices::pdf(file.path(results_dir, "BNST_meta_analysis_forest_plot_recomputed.pdf"), width = 6, height = 5)
metafor::forest(res, slab = meta_input$study, xlab = "Standardized mean difference", mlab = "Random-effects model")
grDevices::dev.off()

message("Wrote prior table and meta-analysis summary to: ", results_dir)
