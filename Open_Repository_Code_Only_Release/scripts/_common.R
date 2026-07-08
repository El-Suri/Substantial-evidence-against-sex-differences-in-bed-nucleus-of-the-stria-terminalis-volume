find_repo_root <- function(start_dir = getwd()) {
  path <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DATA_USE.md")) && dir.exists(file.path(path, "data"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root containing DATA_USE.md and data/.", call. = FALSE)
    }
    path <- parent
  }
}

script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

read_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_public_data <- function(root) {
  public_path <- file.path(root, "data", "open_derived", "bnst_open_derived.csv")
  if (!file.exists(public_path)) {
    stop(
      "Missing open-derived data table: ", public_path, "\n",
      "The code-only release intentionally does not include this file. Add it locally after HCP-approved study IDs are available.",
      call. = FALSE
    )
  }
  read.csv(public_path, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_model_data <- function(public_data, restricted_data) {
  required_public <- c(
    "study_id", "Gender", "TBV_mm3",
    "BNST_left_mm3", "BNST_right_mm3", "BNST_bilateral_mm3"
  )
  missing_public <- setdiff(required_public, names(public_data))
  if (length(missing_public) > 0) {
    stop("Public table is missing: ", paste(missing_public, collapse = ", "), call. = FALSE)
  }

  required_restricted <- c("study_id", "Age_in_Yrs", "Family_ID")
  missing_restricted <- setdiff(required_restricted, names(restricted_data))
  if (length(missing_restricted) > 0) {
    stop("Restricted covariate file is missing: ", paste(missing_restricted, collapse = ", "), call. = FALSE)
  }

  if (!("monozygotic" %in% names(restricted_data))) {
    if ("ZygositySR" %in% names(restricted_data)) {
      restricted_data$monozygotic <- ifelse(restricted_data$ZygositySR == "MZ", 1, 0)
    } else {
      stop("Restricted covariate file must include monozygotic or ZygositySR.", call. = FALSE)
    }
  }

  public_data$study_id <- as.character(public_data$study_id)
  restricted_data$study_id <- as.character(restricted_data$study_id)

  model_data <- merge(public_data, restricted_data, by = "study_id", all = FALSE)
  if (nrow(model_data) != nrow(public_data)) {
    stop(
      "Restricted covariate file did not match all open-derived rows. Matched ",
      nrow(model_data), " of ", nrow(public_data), " rows.",
      call. = FALSE
    )
  }

  model_data$Age_in_Yrs <- as.numeric(model_data$Age_in_Yrs)
  model_data$TBV_mm3 <- as.numeric(model_data$TBV_mm3)
  model_data$BNST_left_mm3 <- as.numeric(model_data$BNST_left_mm3)
  model_data$BNST_right_mm3 <- as.numeric(model_data$BNST_right_mm3)
  model_data$BNST_bilateral_mm3 <- as.numeric(model_data$BNST_bilateral_mm3)
  model_data$monozygotic <- as.integer(model_data$monozygotic)
  model_data$Family_ID <- as.factor(model_data$Family_ID)

  model_data$Age_in_Yrs_Z <- as.numeric(scale(model_data$Age_in_Yrs))
  model_data$TBV_Z <- as.numeric(scale(model_data$TBV_mm3))
  model_data$BNST_bilateral_Z <- as.numeric(scale(model_data$BNST_bilateral_mm3))
  model_data$BNST_left_Z <- as.numeric(scale(model_data$BNST_left_mm3))
  model_data$BNST_right_Z <- as.numeric(scale(model_data$BNST_right_mm3))
  model_data$Gender_dummy <- ifelse(model_data$Gender == "M", 0.5, -0.5)

  model_data <- model_data[order(model_data$Family_ID, model_data$study_id), ]
  mz_pair_id <- character(nrow(model_data))
  for (family in unique(model_data$Family_ID)) {
    idx <- which(model_data$Family_ID == family)
    rows <- model_data[idx, , drop = FALSE]
    mz_pair_id[idx] <- ifelse(
      rows$monozygotic == 1,
      paste0(rows$Family_ID, "_MZpair"),
      paste0(rows$Family_ID, "_", seq_along(idx), "_indep")
    )
  }
  model_data$MZ_pair_ID <- mz_pair_id
  model_data$MZ_pair_ID <- as.factor(model_data$MZ_pair_ID)

  model_data
}

savage_dickey_bf01 <- function(draws, prior_mean, prior_sd) {
  posterior_density <- polspline::logspline(draws)
  polspline::dlogspline(0, posterior_density) / stats::dnorm(0, prior_mean, prior_sd)
}
