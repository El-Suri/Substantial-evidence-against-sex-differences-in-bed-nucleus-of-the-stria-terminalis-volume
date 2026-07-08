args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
this_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/99_release_checks.R"
source(file.path(dirname(normalizePath(this_file)), "_common.R"))

root <- find_repo_root(script_dir())
all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
relative <- sub(paste0("^", normalizePath(root), "/?"), "", normalizePath(all_files, mustWork = FALSE))
is_file <- file.info(all_files)$isdir == FALSE

failures <- character()

bad_ext <- grepl("\\.(rds|RDS|RData|Rhistory)$", relative)
if (any(bad_ext)) failures <- c(failures, paste("Forbidden saved-object/history files:", paste(relative[bad_ext], collapse = ", ")))

image_like <- grepl("\\.nii(\\.gz)?$", relative, ignore.case = TRUE) |
  grepl("(^|/)(T1w|Practice_data|Original|Assessment|restricted_fast_outputs)(/|$)", relative) |
  grepl("T1w_acpc|restore_brain|pve_[0-9]", basename(relative), ignore.case = TRUE)
if (any(image_like)) failures <- c(failures, paste("Image or FAST-derived files are not allowed in the code-only release:", paste(relative[image_like], collapse = ", ")))

forbidden_data <- (
  grepl("^data/open_derived/.*[.]csv$", relative) |
    grepl("^data/segmentations/.*[.]nii([.]gz)?$", relative, ignore.case = TRUE) |
    grepl("^data/prior_sources/.*[.](csv|tsv|xlsx|rds|RDS)$", relative, ignore.case = TRUE) |
    grepl("^results/.*[.](csv|tsv|pdf|rds|RDS|RData)$", relative, ignore.case = TRUE) |
    (grepl("^data/local_inputs/", relative) & !grepl("^data/local_inputs/README[.]md$", relative))
) & is_file
if (any(forbidden_data)) {
  failures <- c(failures, paste("Participant-level/local/result data are not allowed in this code-only release:", paste(relative[forbidden_data], collapse = ", ")))
}

allowed_csv <- "data/restricted_template/restricted_covariates_template.csv"
csv_files <- relative[is_file & grepl("[.]csv$", relative, ignore.case = TRUE)]
unexpected_csv <- setdiff(csv_files, allowed_csv)
if (length(unexpected_csv) > 0) {
  failures <- c(failures, paste("Unexpected CSV files in code-only release:", paste(unexpected_csv, collapse = ", ")))
}

text_ext <- grepl("\\.(R|Rmd|md|txt|csv|tsv|json|yml|yaml|gitignore)$", relative, ignore.case = TRUE)
text_files <- all_files[text_ext & is_file]
for (file in text_files) {
  lines <- tryCatch(readLines(file, warn = FALSE), error = function(e) character())
  local_patterns <- c(
    paste0("/", "Users", "/", "samuelberry"),
    paste0("Cloud", "Storage"),
    paste0("One", "Drive"),
    paste("Royal", "Holloway", "Dropbox")
  )
  hit <- grepl(paste(local_patterns, collapse = "|"), lines)
  if (any(hit)) {
    failures <- c(failures, paste("Local absolute path found in", sub(paste0("^", root, "/?"), "", file)))
  }
}

template <- file.path(root, allowed_csv)
if (!file.exists(template)) {
  failures <- c(failures, paste("Missing restricted covariate template:", allowed_csv))
} else if (length(readLines(template, warn = FALSE)) != 1) {
  failures <- c(failures, "Restricted covariate template must contain only the header row.")
}

if (length(failures) > 0) {
  stop(paste(c("Code-only release checks failed:", failures), collapse = "\n- "), call. = FALSE)
}

message("Code-only release checks passed.")
message("Files scanned: ", length(all_files))
