args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/00_calculate_tbv_from_fast_pves.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
require_packages(c("RNifti"))

subjects_path <- read_arg("subjects", file.path(root, "data", "local_inputs", "study_id_mapping.csv"))
fast_root <- read_arg("fast_root", file.path(root, "data", "restricted_fast_outputs"))
out_path <- read_arg("out", file.path(root, "results", "tbv_from_fast_recomputed.csv"))
ensure_dir(dirname(out_path))

if (!file.exists(subjects_path)) stop("Missing subject list: ", subjects_path, call. = FALSE)
if (!dir.exists(fast_root)) {
  stop(
    "FAST partial-volume output directory not found: ", fast_root, "\n",
    "Provide --fast_root=/path/to/local/FAST_outputs. This directory is not included in the public release.",
    call. = FALSE
  )
}

subjects <- read.csv(subjects_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"study_id" %in% names(subjects)) stop("Mapping file must contain a study_id column.", call. = FALSE)
if (!"hcp_subject_id" %in% names(subjects)) {
  message("No hcp_subject_id column found; using study_id as the local FAST file lookup key.")
  subjects$hcp_subject_id <- subjects$study_id
}
study_ids <- as.character(subjects$study_id)
lookup_ids <- as.character(subjects$hcp_subject_id)

gm_files <- list.files(fast_root, pattern = "pve_1[.]nii([.]gz)?$", recursive = TRUE, full.names = TRUE)
wm_files <- list.files(fast_root, pattern = "pve_2[.]nii([.]gz)?$", recursive = TRUE, full.names = TRUE)
if (length(gm_files) == 0 || length(wm_files) == 0) {
  stop("Could not find FAST *_pve_1.nii(.gz) and *_pve_2.nii(.gz) files in ", fast_root, call. = FALSE)
}

find_subject_file <- function(subject, files, tissue_name) {
  hits <- files[grepl(subject, files, fixed = TRUE)]
  if (length(hits) != 1) {
    stop(
      "Expected exactly one ", tissue_name, " partial-volume file for subject ", subject,
      ", found ", length(hits), ".",
      call. = FALSE
    )
  }
  hits
}

tissue_volume_mm3 <- function(file) {
  img <- RNifti::readNifti(file)
  voxel_dims <- RNifti::pixdim(img)
  voxel_volume <- prod(voxel_dims[seq_len(min(3, length(voxel_dims)))])
  sum(as.numeric(img), na.rm = TRUE) * voxel_volume
}

tbv <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  gm_file <- find_subject_file(lookup_ids[i], gm_files, "grey-matter")
  wm_file <- find_subject_file(lookup_ids[i], wm_files, "white-matter")
  gm_volume <- tissue_volume_mm3(gm_file)
  wm_volume <- tissue_volume_mm3(wm_file)
  data.frame(
    study_id = study_ids[i],
    GM_volume_mm3 = gm_volume,
    WM_volume_mm3 = wm_volume,
    TBV_mm3 = gm_volume + wm_volume,
    stringsAsFactors = FALSE
  )
}))

write.csv(tbv, out_path, row.names = FALSE)
message("Wrote FAST-derived TBV table to: ", out_path)
message("Rows: ", nrow(tbv))
