args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/05_unilateral_sensitivity_and_descriptives.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("brms", "posterior", "polspline", "bayesplot", "performance"))

restricted_path <- read_arg("restricted")
if (is.null(restricted_path)) {
  stop("Provide HCP restricted covariates with --restricted=path/to/restricted_covariates.csv.", call. = FALSE)
}
if (!file.exists(restricted_path)) stop("Restricted covariate file not found: ", restricted_path, call. = FALSE)

iter <- as.integer(read_arg("iter", "10000"))
chains <- as.integer(read_arg("chains", "4"))
cores <- as.integer(read_arg("cores", as.character(chains)))
seed <- as.integer(read_arg("seed", "1234"))

public <- read_public_data(root)
restricted <- read.csv(restricted_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- prepare_model_data(public, restricted)

fmt_mean_sd <- function(x, digits = 2) {
  paste0(round(mean(x, na.rm = TRUE), digits), " (", round(stats::sd(x, na.rm = TRUE), digits), ")")
}

fmt_mean_sd_range <- function(x, digits = 2) {
  paste0(
    round(mean(x, na.rm = TRUE), digits),
    " (", round(stats::sd(x, na.rm = TRUE), digits), ") [",
    round(min(x, na.rm = TRUE), digits), "-",
    round(max(x, na.rm = TRUE), digits), "]"
  )
}

classify_family_type <- function(data) {
  if (!("ZygositySR" %in% names(data))) {
    return(rep("not classified", nrow(data)))
  }
  family_type <- character(nrow(data))
  for (family in unique(data$Family_ID)) {
    idx <- which(data$Family_ID == family)
    z <- data$ZygositySR[idx]
    label <- if (length(idx) == 2 && sum(z == "MZ", na.rm = TRUE) == 2) {
      "MZ twin pair"
    } else if (length(idx) == 2 && sum(z == "NotMZ", na.rm = TRUE) == 2) {
      "DZ twin pair"
    } else if (length(idx) >= 2 && all(z == "NotTwin", na.rm = TRUE)) {
      "non-twin siblings"
    } else {
      "unrelated individual"
    }
    family_type[idx] <- label
  }
  family_type
}

df$Sex <- ifelse(df$Gender == "M", "Male", "Female")
df$family_type <- classify_family_type(df)

descriptives_by_sex <- lapply(c("Total", "Female", "Male"), function(group) {
  subset <- if (group == "Total") df else df[df$Sex == group, , drop = FALSE]
  data.frame(
    Sex = group,
    N = nrow(subset),
    age = fmt_mean_sd_range(subset$Age_in_Yrs),
    TBV = fmt_mean_sd(subset$TBV_mm3, 1),
    bilateral_BNST = fmt_mean_sd(subset$BNST_bilateral_mm3, 2),
    left_BNST = fmt_mean_sd(subset$BNST_left_mm3, 2),
    right_BNST = fmt_mean_sd(subset$BNST_right_mm3, 2),
    MZ = sum(subset$family_type == "MZ twin pair"),
    DZ = sum(subset$family_type == "DZ twin pair"),
    non_twin_siblings = sum(subset$family_type == "non-twin siblings"),
    unrelated = sum(subset$family_type == "unrelated individual"),
    stringsAsFactors = FALSE
  )
})
descriptives <- do.call(rbind, descriptives_by_sex)

pooled_d <- 0.6591
pooled_se <- 0.2907
prior_sd <- pooled_se * 2
model_prior <- brms::set_prior(
  paste0("normal(", pooled_d, ", ", prior_sd, ")"),
  class = "b",
  coef = "Gender_dummy"
) +
  brms::set_prior("normal(0,1)", class = "b") +
  brms::set_prior("normal(0,1)", class = "sigma")

fit_unilateral <- function(outcome, label) {
  formula <- stats::as.formula(
    paste0(
      outcome,
      " ~ Gender_dummy + Age_in_Yrs_Z + I(Age_in_Yrs_Z^2) + TBV_Z + ",
      "(1 | Family_ID) + (1 | MZ_pair_ID)"
    )
  )
  fit <- brms::brm(
    formula = formula,
    data = df,
    prior = model_prior,
    sample_prior = "yes",
    control = list(max_treedepth = 12, adapt_delta = 0.99),
    iter = iter,
    cores = cores,
    chains = chains,
    seed = seed,
    family = gaussian(),
    refresh = 0
  )
  fixed <- as.data.frame(brms::fixef(fit, probs = c(0.025, 0.975)))
  draws <- posterior::as_draws_df(fit)
  r2_bayes <- performance::r2(fit)
  diagnostics <- as.data.frame(posterior::summarise_draws(fit))
  sampler <- bayesplot::nuts_params(fit)

  data.frame(
    hemisphere = label,
    sex_beta = fixed["Gender_dummy", "Estimate"],
    sex_ci_low = fixed["Gender_dummy", "Q2.5"],
    sex_ci_high = fixed["Gender_dummy", "Q97.5"],
    BF01 = savage_dickey_bf01(draws$b_Gender_dummy, pooled_d, prior_sd),
    TBV_beta = fixed["TBV_Z", "Estimate"],
    TBV_ci_low = fixed["TBV_Z", "Q2.5"],
    TBV_ci_high = fixed["TBV_Z", "Q97.5"],
    age_beta = fixed["Age_in_Yrs_Z", "Estimate"],
    age_ci_low = fixed["Age_in_Yrs_Z", "Q2.5"],
    age_ci_high = fixed["Age_in_Yrs_Z", "Q97.5"],
    marginal_R2 = unname(r2_bayes$R2_Bayes_marginal[1]),
    conditional_R2 = unname(r2_bayes$R2_Bayes[1]),
    max_rhat = max(diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
    divergent_transitions = sum(sampler$Parameter == "divergent__" & sampler$Value == 1, na.rm = TRUE)
  )
}

unilateral <- rbind(
  fit_unilateral("BNST_right_Z", "Right"),
  fit_unilateral("BNST_left_Z", "Left")
)

paired_lr <- stats::t.test(df$BNST_left_mm3, df$BNST_right_mm3, paired = TRUE)
paired_row <- data.frame(
  hemisphere = "Left minus right",
  sex_beta = NA_real_,
  sex_ci_low = NA_real_,
  sex_ci_high = NA_real_,
  BF01 = NA_real_,
  TBV_beta = NA_real_,
  TBV_ci_low = NA_real_,
  TBV_ci_high = NA_real_,
  age_beta = NA_real_,
  age_ci_low = NA_real_,
  age_ci_high = NA_real_,
  marginal_R2 = NA_real_,
  conditional_R2 = NA_real_,
  max_rhat = NA_real_,
  min_bulk_ess = NA_real_,
  min_tail_ess = NA_real_,
  divergent_transitions = NA_real_,
  paired_t = unname(paired_lr$statistic),
  paired_df = unname(paired_lr$parameter),
  paired_p = paired_lr$p.value,
  paired_mean_difference = unname(paired_lr$estimate)
)

unilateral$paired_t <- NA_real_
unilateral$paired_df <- NA_real_
unilateral$paired_p <- NA_real_
unilateral$paired_mean_difference <- NA_real_
results <- rbind(unilateral, paired_row)

results_dir <- file.path(root, "results")
ensure_dir(results_dir)
write.csv(descriptives, file.path(results_dir, "descriptives_recomputed.csv"), row.names = FALSE)
write.csv(results, file.path(results_dir, "unilateral_results_recomputed.csv"), row.names = FALSE)
message("Wrote unilateral sensitivity and descriptives to: ", results_dir)
