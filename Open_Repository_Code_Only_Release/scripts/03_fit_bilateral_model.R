args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/03_fit_bilateral_model.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("brms", "posterior", "polspline", "bayesplot"))

restricted_path <- read_arg("restricted")
if (is.null(restricted_path)) {
  stop(
    "Provide HCP restricted covariates with --restricted=path/to/restricted_covariates.csv. ",
    "Use data/restricted_template/restricted_covariates_template.csv as the column template.",
    call. = FALSE
  )
}
if (!file.exists(restricted_path)) stop("Restricted covariate file not found: ", restricted_path, call. = FALSE)

iter <- as.integer(read_arg("iter", "10000"))
chains <- as.integer(read_arg("chains", "4"))
cores <- as.integer(read_arg("cores", as.character(chains)))
seed <- as.integer(read_arg("seed", "1234"))

public <- read_public_data(root)
restricted <- read.csv(restricted_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- prepare_model_data(public, restricted)

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

fit <- brms::brm(
  formula = BNST_bilateral_Z ~ Gender_dummy + Age_in_Yrs_Z + I(Age_in_Yrs_Z^2) + TBV_Z +
    (1 | Family_ID) + (1 | MZ_pair_ID),
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
bf01 <- savage_dickey_bf01(draws$b_Gender_dummy, pooled_d, prior_sd)
r2 <- brms::bayes_R2(fit, summary = TRUE)
diagnostics <- as.data.frame(posterior::summarise_draws(fit))
sampler <- bayesplot::nuts_params(fit)

results_dir <- file.path(root, "results")
ensure_dir(results_dir)

main_results <- data.frame(
  term = rownames(fixed),
  estimate = fixed$Estimate,
  est_error = fixed$Est.Error,
  ci_low = fixed$Q2.5,
  ci_high = fixed$Q97.5,
  stringsAsFactors = FALSE
)
main_results$BF01_sex <- NA_real_
main_results$BF01_sex[main_results$term == "Gender_dummy"] <- bf01

write.csv(main_results, file.path(results_dir, "bilateral_model_results.csv"), row.names = FALSE)

diagnostic_results <- data.frame(
  n_subjects = nrow(df),
  prior_mean_sex = pooled_d,
  prior_sd_sex = prior_sd,
  bayes_R2_estimate = r2[1, "Estimate"],
  bayes_R2_error = r2[1, "Est.Error"],
  max_rhat = max(diagnostics$rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
  divergent_transitions = sum(sampler$Parameter == "divergent__" & sampler$Value == 1, na.rm = TRUE)
)
write.csv(diagnostic_results, file.path(results_dir, "bilateral_model_diagnostics.csv"), row.names = FALSE)

message("Wrote bilateral model results to: ", results_dir)
message("BF01 for sex effect: ", signif(bf01, 4))
