args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/03_fit_bilateral_model.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c(
  "brms", "posterior", "polspline", "bayesplot", "ggplot2",
  "dplyr", "tidyr", "tidybayes", "patchwork", "loo"
))

restricted_path <- read_arg("restricted")
if (is.null(restricted_path)) {
  stop(
    "Provide HCP restricted covariates with --restricted=data/local_inputs/restricted_covariates.csv. ",
    "Use data/restricted_template/restricted_covariates_template.csv as the column template.",
    call. = FALSE
  )
}
if (!file.exists(restricted_path)) stop("Restricted covariate file not found: ", restricted_path, call. = FALSE)

iter <- as.integer(read_arg("iter", "10000"))
chains <- as.integer(read_arg("chains", "4"))
cores <- as.integer(read_arg("cores", as.character(chains)))
seed <- as.integer(read_arg("seed", "1234"))

results_dir <- file.path(root, "results")
figures_dir <- file.path(root, "figures")
ensure_dir(results_dir)
ensure_dir(figures_dir)

# -------------------------------------------------------------------------
# Load and prepare the analysis data.
# -------------------------------------------------------------------------
# The public table uses study_id, not HCP Subject ID. Exact age and family
# structure are read from a local restricted covariate file that is never
# uploaded. prepare_model_data() standardizes age, TBV, and BNST volumes and
# creates the MZ-pair random-effect grouping variable used in the manuscript.

public <- read_public_data(root)
restricted <- read.csv(restricted_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- prepare_model_data(public, restricted)

# Literature-informed prior from scripts/02_create_prior_table.R.
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

# -------------------------------------------------------------------------
# Main manuscript model: bilateral BNST volume adjusted for age, nonlinear age,
# TBV, family, and MZ-pair structure.
# -------------------------------------------------------------------------

main_formula <- BNST_bilateral_Z ~ Gender_dummy + Age_in_Yrs_Z + I(Age_in_Yrs_Z^2) + TBV_Z +
  (1 | Family_ID) + (1 | MZ_pair_ID)

bayesian_model <- brms::brm(
  formula = main_formula,
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

# The no-TBV model is not the primary inference, but it documents the model
# comparison and makes clear how TBV changes the sex coefficient.
notbv_formula <- BNST_bilateral_Z ~ Gender_dummy + Age_in_Yrs_Z + I(Age_in_Yrs_Z^2) +
  (1 | Family_ID) + (1 | MZ_pair_ID)

bayesian_model_notbv <- brms::brm(
  formula = notbv_formula,
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

# -------------------------------------------------------------------------
# Numeric summaries and diagnostics.
# -------------------------------------------------------------------------

fixed <- as.data.frame(brms::fixef(bayesian_model, probs = c(0.025, 0.975)))
draws <- posterior::as_draws_df(bayesian_model)
bf01 <- savage_dickey_bf01(draws$b_Gender_dummy, pooled_d, prior_sd)
r2 <- brms::bayes_R2(bayesian_model, summary = TRUE)
diagnostics <- as.data.frame(posterior::summarise_draws(bayesian_model))
sampler <- bayesplot::nuts_params(bayesian_model)

main_results <- data.frame(
  term = rownames(fixed),
  estimate = fixed$Estimate,
  est_error = fixed$Est.Error,
  ci_low = fixed$Q2.5,
  ci_high = fixed$Q97.5,
  BF01_sex = NA_real_,
  stringsAsFactors = FALSE
)
main_results$BF01_sex[main_results$term == "Gender_dummy"] <- bf01
write.csv(main_results, file.path(results_dir, "bilateral_model_results.csv"), row.names = FALSE)

vc <- brms::VarCorr(bayesian_model, summary = TRUE)
sd_family <- vc$Family_ID$sd[1, "Estimate"]
sd_mzpair <- vc$MZ_pair_ID$sd[1, "Estimate"]
sd_resid <- vc$residual__$sd[1, "Estimate"]
var_family <- sd_family^2
var_mzpair <- sd_mzpair^2
var_resid <- sd_resid^2
var_total <- var_family + var_mzpair + var_resid

diagnostic_results <- data.frame(
  n_subjects = nrow(df),
  prior_mean_sex = pooled_d,
  prior_sd_sex = prior_sd,
  BF01_sex = bf01,
  bayes_R2_estimate = r2[1, "Estimate"],
  bayes_R2_error = r2[1, "Est.Error"],
  max_rhat = max(diagnostics$rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
  divergent_transitions = sum(sampler$Parameter == "divergent__" & sampler$Value == 1, na.rm = TRUE),
  family_variance = var_family,
  mz_pair_variance = var_mzpair,
  residual_variance = var_resid,
  icc_family = var_family / var_total,
  icc_family_plus_mz_pair = (var_family + var_mzpair) / var_total,
  stringsAsFactors = FALSE
)
write.csv(diagnostic_results, file.path(results_dir, "bilateral_model_diagnostics.csv"), row.names = FALSE)

loo_tbv <- loo::loo(bayesian_model)
loo_notbv <- loo::loo(bayesian_model_notbv)
loo_comparison <- as.data.frame(loo::loo_compare(loo_tbv, loo_notbv))
loo_comparison$model <- rownames(loo_comparison)
write.csv(loo_comparison, file.path(results_dir, "loo_model_comparison_tbv_vs_no_tbv.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Shared figure styling.
# -------------------------------------------------------------------------

shared_theme <- ggplot2::theme_classic(base_size = 11, base_family = "Helvetica") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    axis.title = ggplot2::element_text(size = 10),
    axis.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

point_col <- "purple"
line_col <- "black"

# -------------------------------------------------------------------------
# Figure: model diagnostics.
# -------------------------------------------------------------------------

fit_mean <- brms::fitted(bayesian_model, summary = TRUE)[, "Estimate"]
res_mean <- brms::residuals(bayesian_model, summary = TRUE, type = "ordinary")[, "Estimate"]
df_resid <- data.frame(fitted = fit_mean, residuals = res_mean)

p_diag1 <- ggplot2::ggplot(df_resid, ggplot2::aes(x = fitted, y = residuals)) +
  ggplot2::geom_point(colour = point_col, alpha = 0.55, size = 1.8) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = line_col, linewidth = 0.7) +
  ggplot2::labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  shared_theme

qq <- stats::qqnorm(res_mean, plot.it = FALSE)
probs <- c(0.25, 0.75)
y_q <- stats::quantile(res_mean, probs, na.rm = TRUE)
x_q <- stats::qnorm(probs)
slope <- diff(y_q) / diff(x_q)
interc <- y_q[1] - slope * x_q[1]
df_qq <- data.frame(theoretical = qq$x, sample = qq$y)

p_diag2 <- ggplot2::ggplot(df_qq, ggplot2::aes(x = theoretical, y = sample)) +
  ggplot2::geom_point(colour = point_col, alpha = 0.55, size = 1.8) +
  ggplot2::geom_abline(intercept = interc, slope = slope, colour = line_col, linetype = "dashed", linewidth = 0.7) +
  ggplot2::labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Sample quantiles") +
  shared_theme

p_diag3 <- brms::pp_check(bayesian_model, ndraws = 100, type = "dens_overlay") +
  ggplot2::labs(title = "Posterior Predictive Check") +
  shared_theme +
  ggplot2::theme(legend.position = "bottom")

combined_diag <- p_diag1 + p_diag2 + p_diag3 +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(
    title = "Bayesian Model Diagnostics",
    subtitle = "Bilateral BNST volume ~ Sex + Age + TBV + (1 | Family_ID) + (1 | MZ_pair_ID)"
  )

ggplot2::ggsave(file.path(figures_dir, "model_diagnostics.pdf"), combined_diag, width = 10, height = 4, device = "pdf")
ggplot2::ggsave(file.path(figures_dir, "model_diagnostics.png"), combined_diag, width = 10, height = 4, dpi = 300)

# -------------------------------------------------------------------------
# Figure: Savage-Dickey BF and posterior predictions by sex.
# -------------------------------------------------------------------------

savage_dickey_plot <- function(x, x_0 = 0, prior_mean = 0, prior_sd = 1) {
  posterior_density <- polspline::logspline(x)
  posterior_w <- polspline::dlogspline(x_0, posterior_density)
  prior_w <- stats::dnorm(x_0, prior_mean, prior_sd)
  bf <- posterior_w / prior_w

  x_seq <- seq(posterior_density$range[1], posterior_density$range[2] + diff(posterior_density$range) / 3, length.out = 500)
  curves <- rbind(
    data.frame(x = x_seq, y = polspline::dlogspline(x_seq, posterior_density), curve = "Posterior"),
    data.frame(x = x_seq, y = stats::dnorm(x_seq, prior_mean, prior_sd), curve = "Prior")
  )
  points <- data.frame(x = c(x_0, x_0), y = c(posterior_w, prior_w), curve = c("Posterior", "Prior"))
  bf_label <- sprintf("atop(BF[0*1] == %.2f, paste('(evidence for ', H[0], ')'))", bf)

  plot <- ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = data.frame(x = x),
      ggplot2::aes(x = x, y = ggplot2::after_stat(density)),
      fill = "#5C6D70", alpha = 0.25, bins = 40, colour = NA
    ) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = x, y = y, colour = curve), linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = x_0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    ggplot2::geom_point(data = points, ggplot2::aes(x = x, y = y, colour = curve), size = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = x_0, xend = x_0, y = posterior_w, yend = prior_w),
                          colour = "grey30", linewidth = 0.5, linetype = "dotted") +
    ggplot2::scale_colour_manual(values = c("Posterior" = "#5C6D70", "Prior" = "#E88873"), name = NULL) +
    ggplot2::annotate("text", x = max(x_seq), y = max(curves$y) * 0.97,
                      label = bf_label, parse = TRUE, hjust = 1, size = 3.2, colour = "grey20") +
    ggplot2::labs(title = "Savage-Dickey Bayes Factor", x = expression(beta[Sex]), y = "Density") +
    shared_theme +
    ggplot2::theme(legend.position = c(0.97, 0.75), legend.justification = c(1, 1))

  list(plot = plot, bf = bf)
}

sd_result <- savage_dickey_plot(draws$b_Gender_dummy, prior_mean = pooled_d, prior_sd = prior_sd)
p_bf1 <- sd_result$plot

ppd <- brms::posterior_predict(bayesian_model)
ppd_df <- as.data.frame(t(ppd))
ppd_df$Gender <- df$Gender_dummy

ppd_long <- ppd_df |>
  tidyr::pivot_longer(cols = -Gender, names_to = "sample", values_to = "Predicted_Value") |>
  dplyr::mutate(
    Predicted_Value = attr(scale(df$BNST_bilateral_mm3), "scaled:center") +
      Predicted_Value * attr(scale(df$BNST_bilateral_mm3), "scaled:scale"),
    plot_col = ifelse(Gender == 0.5, "M", "F")
  )

df_long <- df |>
  dplyr::mutate(sample = dplyr::row_number()) |>
  tidyr::pivot_longer(cols = BNST_bilateral_Z, names_to = "Variable", values_to = "Observed_Value") |>
  dplyr::mutate(
    Observed_Value = attr(scale(df$BNST_bilateral_mm3), "scaled:center") +
      Observed_Value * attr(scale(df$BNST_bilateral_mm3), "scaled:scale"),
    plot_col = ifelse(Gender_dummy == 0.5, "M", "F")
  )

p_bf2 <- ppd_long |>
  ggplot2::ggplot(ggplot2::aes(x = Gender, y = Predicted_Value, colour = plot_col, fill = plot_col)) +
  tidybayes::stat_halfeye(
    mapping = ggplot2::aes(thickness = ggplot2::after_stat(pdf * n)),
    .width = 0.95,
    point_interval = "mean_hdi",
    interval_color = "black",
    point_color = "black",
    justification = 0,
    scale = 0.5,
    linewidth = 8,
    size = 0
  ) +
  ggplot2::geom_point(
    data = df_long,
    ggplot2::aes(y = Observed_Value, x = Gender_dummy - 0.16),
    position = ggplot2::position_jitter(width = 0.1),
    size = 1.5,
    alpha = 0.8
  ) +
  ggplot2::scale_colour_manual(values = c("F" = "#E88873", "M" = "#5C6D70"), guide = "none") +
  ggplot2::scale_fill_manual(values = c("F" = "#E88873", "M" = "#5C6D70"), guide = "none") +
  ggplot2::coord_flip(xlim = c(-1, 1)) +
  ggplot2::scale_x_continuous(breaks = c(-0.5, 0.5), labels = c("F", "M")) +
  ggplot2::labs(title = "Posterior Predictions by Sex", x = "Sex", y = "BNST volume (mm³)") +
  shared_theme

combined_bf <- p_bf1 + p_bf2 +
  patchwork::plot_layout(ncol = 2, widths = c(1.1, 1)) +
  patchwork::plot_annotation(title = "Sex effect on Bilateral BNST Volume")

ggplot2::ggsave(file.path(figures_dir, "BNST_gender_BF_ppd.pdf"), combined_bf, width = 8, height = 4, device = "pdf")
ggplot2::ggsave(file.path(figures_dir, "BNST_gender_BF_ppd.png"), combined_bf, width = 8, height = 4, dpi = 300)

# -------------------------------------------------------------------------
# Figure: raw observed versus TBV-adjusted BNST volume by sex.
# -------------------------------------------------------------------------

bnst_center <- attr(scale(df$BNST_bilateral_mm3), "scaled:center")
bnst_scale <- attr(scale(df$BNST_bilateral_mm3), "scaled:scale")
df$BNST_raw <- df$BNST_bilateral_mm3
df$Sex <- ifelse(df$Gender_dummy == -0.5, "Female", "Male")

beta_TBV <- draws$b_TBV_Z
df$BNST_adj_Z <- df$BNST_bilateral_Z - mean(beta_TBV) * df$TBV_Z
df$BNST_adj_mm3 <- bnst_center + df$BNST_adj_Z * bnst_scale

sex_colours <- c("Female" = "#E88873", "Male" = "#5C6D70")

p_tbv1 <- ggplot2::ggplot(df, ggplot2::aes(x = Sex, y = BNST_raw, colour = Sex, fill = Sex)) +
  ggplot2::geom_violin(alpha = 0.25, trim = FALSE, linewidth = 0.6) +
  ggplot2::geom_jitter(width = 0.08, height = 0, size = 1.5, alpha = 0.7) +
  ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.25, linewidth = 0.6, colour = "grey20", fatten = 1) +
  ggplot2::scale_colour_manual(values = sex_colours, guide = "none") +
  ggplot2::scale_fill_manual(values = sex_colours, guide = "none") +
  ggplot2::labs(title = "Observed", x = NULL, y = expression(BNST~volume~(mm^3))) +
  shared_theme

p_tbv2 <- ggplot2::ggplot(df, ggplot2::aes(x = Sex, y = BNST_adj_mm3, colour = Sex, fill = Sex)) +
  ggplot2::geom_violin(alpha = 0.25, trim = FALSE, linewidth = 0.6) +
  ggplot2::geom_jitter(width = 0.08, height = 0, size = 1.5, alpha = 0.7) +
  ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.25, linewidth = 0.6, colour = "grey20", fatten = 1) +
  ggplot2::scale_colour_manual(values = sex_colours, guide = "none") +
  ggplot2::scale_fill_manual(values = sex_colours, guide = "none") +
  ggplot2::labs(title = "TBV-adjusted", x = NULL, y = expression(BNST~volume~(mm^3)~"(TBV-adjusted)")) +
  shared_theme

combined_tbv <- p_tbv1 + p_tbv2 +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(title = "Bilateral BNST Volume by Sex")

ggplot2::ggsave(file.path(figures_dir, "BNST_raw_vs_adjusted.pdf"), combined_tbv, width = 8, height = 4, device = "pdf")
ggplot2::ggsave(file.path(figures_dir, "BNST_raw_vs_adjusted.png"), combined_tbv, width = 8, height = 4, dpi = 300)

message("Wrote bilateral model results to: ", results_dir)
message("Wrote manuscript figures to: ", figures_dir)
message("BF01 for sex effect: ", signif(bf01, 4))
