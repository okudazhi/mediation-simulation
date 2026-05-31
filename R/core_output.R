################################################################################
# core_output.R
#
# Lightweight file-system helpers shared by the simulation and real-data paths.
################################################################################

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

write_csv_safe <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(x, file = path, row.names = row.names, na = "")
  invisible(path)
}

save_rds_safe <- function(x, path) {
  ensure_dir(dirname(path))
  saveRDS(x, file = path)
  invisible(path)
}

write_lines_safe <- function(lines, path) {
  ensure_dir(dirname(path))
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

new_run_dir <- function(base_dir, prefix = "run") {
  ensure_dir(base_dir)
  out <- file.path(base_dir, paste0(prefix, "_", project_timestamp()))
  ensure_dir(out)
  out
}
