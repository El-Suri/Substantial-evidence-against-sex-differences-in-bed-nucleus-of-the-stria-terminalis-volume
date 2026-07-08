args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/01_extract_bnst_volumes.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("RNifti"))

seg_dir <- read_arg("seg_dir", file.path(root, "data", "segmentations"))
out_path <- read_arg("out", file.path(root, "results", "extracted_bnst_volumes.csv"))
ensure_dir(dirname(out_path))

files <- list.files(seg_dir, pattern = "\\.nii\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("No .nii.gz masks found in ", seg_dir, call. = FALSE)

results <- lapply(files, function(file) {
  img <- RNifti::readNifti(file)
  voxel_dims <- RNifti::pixdim(img)
  voxel_size <- prod(voxel_dims[seq_len(min(3, length(voxel_dims)))])
  right <- sum(img == 1) * voxel_size
  left <- sum(img == 2) * voxel_size
  data.frame(
    study_id = sub("\\.nii\\.gz$", "", basename(file)),
    BNST_left_mm3 = left,
    BNST_right_mm3 = right,
    BNST_bilateral_mm3 = left + right,
    stringsAsFactors = FALSE
  )
})

volumes <- do.call(rbind, results)
volumes <- volumes[order(volumes$study_id), ]
write.csv(volumes, out_path, row.names = FALSE)

public_path <- file.path(root, "data", "open_derived", "bnst_open_derived.csv")
if (file.exists(public_path)) {
  public <- read_public_data(root)
  merged <- merge(public, volumes, by = "study_id", suffixes = c("_public", "_extracted"))
  max_diff <- max(abs(as.numeric(merged$BNST_bilateral_mm3_public) - merged$BNST_bilateral_mm3_extracted))
  message("Maximum bilateral-volume difference versus open-derived table: ", signif(max_diff, 4), " mm3")
} else {
  message("No open-derived table found; skipped volume comparison.")
}

message("Wrote: ", out_path)
message("Extracted masks: ", nrow(volumes))
