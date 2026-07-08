args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/00_prepare_open_derived_data.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("RNifti"))

demographics_path <- read_arg("demographics", file.path(root, "data", "local_inputs", "open_demographics.csv"))
tbv_path <- read_arg("tbv", file.path(root, "data", "local_inputs", "tbv_source.csv"))
seg_dir <- read_arg("seg_dir", file.path(root, "data", "segmentations"))
out_path <- read_arg("out", file.path(root, "results", "bnst_open_derived_recomputed.csv"))
ensure_dir(dirname(out_path))

for (path in c(demographics_path, tbv_path, seg_dir)) {
  if (!file.exists(path) && !dir.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}

demographics <- read.csv(demographics_path, stringsAsFactors = FALSE, check.names = FALSE)
tbv <- read.csv(tbv_path, stringsAsFactors = FALSE, check.names = FALSE)

missing_demo <- setdiff(c("study_id", "Gender", "Age_group"), names(demographics))
if (length(missing_demo) > 0) {
  stop("Demographics file is missing: ", paste(missing_demo, collapse = ", "), call. = FALSE)
}

missing_tbv <- setdiff(c("study_id", "TBV_mm3"), names(tbv))
if (length(missing_tbv) > 0) {
  stop("TBV source file is missing: ", paste(missing_tbv, collapse = ", "), call. = FALSE)
}

files <- list.files(seg_dir, pattern = "[.]nii[.]gz$", full.names = TRUE)
if (length(files) == 0) stop("No .nii.gz masks found in ", seg_dir, call. = FALSE)

volumes <- do.call(rbind, lapply(files, function(file) {
  img <- RNifti::readNifti(file)
  voxel_dims <- RNifti::pixdim(img)
  voxel_size <- prod(voxel_dims[seq_len(min(3, length(voxel_dims)))])
  right <- sum(img == 1) * voxel_size
  left <- sum(img == 2) * voxel_size
  data.frame(
    study_id = sub("[.]nii[.]gz$", "", basename(file)),
    BNST_left_mm3 = left,
    BNST_right_mm3 = right,
    BNST_bilateral_mm3 = left + right,
    segmentation_file = file.path("data", "segmentations", basename(file)),
    stringsAsFactors = FALSE
  )
}))

open_data <- merge(demographics, tbv, by = "study_id", all = FALSE)
open_data <- merge(open_data, volumes, by = "study_id", all = FALSE)

expected_n <- nrow(demographics)
if (nrow(open_data) != expected_n) {
  stop(
    "Only matched ", nrow(open_data), " of ", expected_n,
    " demographic rows across demographics, TBV, and masks.",
    call. = FALSE
  )
}

open_data$TBV_mm3 <- as.numeric(open_data$TBV_mm3)
open_data$BNST_left_mm3 <- as.numeric(open_data$BNST_left_mm3)
open_data$BNST_right_mm3 <- as.numeric(open_data$BNST_right_mm3)
open_data$BNST_bilateral_mm3 <- as.numeric(open_data$BNST_bilateral_mm3)

open_data$BNST_left_div_TBV <- open_data$BNST_left_mm3 / open_data$TBV_mm3
open_data$BNST_right_div_TBV <- open_data$BNST_right_mm3 / open_data$TBV_mm3
open_data$BNST_bilateral_div_TBV <- open_data$BNST_bilateral_mm3 / open_data$TBV_mm3

open_data <- open_data[order(open_data$study_id), c(
  "study_id", "Gender", "Age_group", "TBV_mm3",
  "BNST_left_mm3", "BNST_right_mm3", "BNST_bilateral_mm3",
  "BNST_left_div_TBV", "BNST_right_div_TBV", "BNST_bilateral_div_TBV",
  "segmentation_file"
)]

write.csv(open_data, out_path, row.names = FALSE)
message("Wrote recomputed open-derived table to: ", out_path)
message("Rows: ", nrow(open_data))
