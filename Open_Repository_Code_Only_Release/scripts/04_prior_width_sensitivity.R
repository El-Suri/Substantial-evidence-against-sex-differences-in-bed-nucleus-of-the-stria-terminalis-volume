args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/04_prior_width_sensitivity.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("brms", "posterior", "polspline"))

restricted_path <- read_arg("restricted")
if (is.null(restricted_path)) {
  stop("Provide HCP restricted covariates with --restricted=path/to/restricted_covariates.csv.", call. = FALSE)
}
if (!file.exists(restricted_path)) stop("Restricted covariate file not found: ", restricted_path, call. = FALSE)

iter <- as.integer(read_arg("iter", "10000"))
chains <- as.integer(read_arg("chains", "4"))
cores <- as.integer(read_arg("cores", as.character(chains)))
seed <- as.integer(read_arg("seed", "1234"))
n_widths <- as.integer(read_arg("n_widths", "25"))

public <- read_public_data(root)
restricted <- read.csv(restricted_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- prepare_model_data(public, restricted)

pooled_d <- 0.6591
pooled_se <- 0.2907
prior_widths <- seq(1, 8, length.out = n_widths)

fit_one_width <- function(width) {
  prior_sd <- pooled_se * width
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
  data.frame(
    prior_SD_scaling = width,
    prior_sd = prior_sd,
    Estimate = fixed["Gender_dummy", "Estimate"],
    Est.Error = fixed["Gender_dummy", "Est.Error"],
    CI.Lower = fixed["Gender_dummy", "Q2.5"],
    CI.Upper = fixed["Gender_dummy", "Q97.5"],
    BF01 = savage_dickey_bf01(draws$b_Gender_dummy, pooled_d, prior_sd)
  )
}

results <- do.call(rbind, lapply(prior_widths, fit_one_width))
results_dir <- file.path(root, "results")
ensure_dir(results_dir)
write.csv(results, file.path(results_dir, "BF_bilateral_prior_width_sensitivity_recomputed.csv"), row.names = FALSE)
message("Wrote prior-width sensitivity results to: ", results_dir)

