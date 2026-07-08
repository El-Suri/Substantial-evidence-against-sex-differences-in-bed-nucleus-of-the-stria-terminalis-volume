args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/04_prior_width_sensitivity.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("brms", "posterior", "polspline"))

restricted_path <- read_arg("restricted")
if (is.null(restricted_path)) {
  stop("Provide HCP restricted covariates with --restricted=data/local_inputs/restricted_covariates.csv.", call. = FALSE)
}
if (!file.exists(restricted_path)) stop("Restricted covariate file not found: ", restricted_path, call. = FALSE)

iter <- as.integer(read_arg("iter", "10000"))
chains <- as.integer(read_arg("chains", "4"))
cores <- as.integer(read_arg("cores", as.character(chains)))
seed <- as.integer(read_arg("seed", "1234"))
n_widths <- as.integer(read_arg("n_widths", "25"))

results_dir <- file.path(root, "results")
figures_dir <- file.path(root, "figures")
ensure_dir(results_dir)
ensure_dir(figures_dir)

public <- read_public_data(root)
restricted <- read.csv(restricted_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- prepare_model_data(public, restricted)

# The manuscript tests whether the conclusion is stable as the prior is made
# progressively wider. Each point below is a full refit of the same Bayesian
# mixed-effects model with a different SD on the sex-effect prior.
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
write.csv(results, file.path(results_dir, "BF_bilateral_prior_width_sensitivity.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Recreate the manuscript prior-width sensitivity figure.
# -------------------------------------------------------------------------

line_error_CI <- function(x, y_err, col, alphabar) {
  rgbbarcol <- grDevices::col2rgb(col) / 255
  graphics::polygon(
    c(x, rev(x)),
    c(y_err[1, ], rev(y_err[2, ])),
    col = grDevices::rgb(rgbbarcol[1], rgbbarcol[2], rgbbarcol[3], alphabar),
    border = NA
  )
}

plot_prior_sensitivity <- function(res_all) {
  bf_ylim <- range(c(1, sqrt(10), res_all$BF01), na.rm = TRUE)
  beta_ylim <- range(c(res_all$CI.Lower, res_all$CI.Upper, 0), na.rm = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), lwd = 1, mar = c(3.2, 3.2, 1.2, 0.8), mgp = c(1.8, 0.5, 0))

  graphics::plot(
    res_all$prior_sd, res_all$BF01,
    type = "o",
    ylim = bf_ylim,
    ylab = expression(BF["01"]),
    pch = 19,
    col = "lightblue",
    xlab = "Prior SD"
  )
  graphics::abline(h = 1, lty = 2)
  graphics::abline(h = sqrt(10), lty = 3)
  graphics::mtext("a", side = 3, line = 0.2, adj = 0, font = 2)

  graphics::plot(
    res_all$prior_sd, res_all$Estimate,
    type = "o",
    ylim = beta_ylim,
    ylab = expression(beta[Sex]),
    pch = 19,
    col = "lightblue",
    xlab = "Prior SD"
  )
  line_error_CI(res_all$prior_sd, rbind(res_all$CI.Lower, res_all$CI.Upper), col = "lightblue", alphabar = 0.25)
  graphics::lines(res_all$prior_sd, res_all$Estimate, type = "o", pch = 19, col = "lightblue")
  graphics::abline(h = 0, lty = 2)
  graphics::mtext("b", side = 3, line = 0.2, adj = 0, font = 2)
}

grDevices::pdf(file.path(figures_dir, "BF_bilateral_prior_width_sensitivity.pdf"), width = 8, height = 4)
plot_prior_sensitivity(results)
grDevices::dev.off()

grDevices::png(file.path(figures_dir, "BF_bilateral_prior_width_sensitivity.png"), width = 2400, height = 1200, res = 300)
plot_prior_sensitivity(results)
grDevices::dev.off()

message("Wrote prior-width sensitivity results to: ", results_dir)
message("Wrote prior-width sensitivity figure to: ", figures_dir)
