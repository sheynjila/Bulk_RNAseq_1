# ==============================================================================
# ORGANIZE PROJECT FILES
# Run 01_create_project_directories.R first, then run this script from the
# project root folder. Existing destination files are never overwritten.
# This organizer stays in the project root so it can be run again; a copy is
# also placed in the scripts directory.
# ==============================================================================

required_directories <- c("scripts", "data", "plots")
missing_directories <- required_directories[!dir.exists(required_directories)]

if (length(missing_directories) > 0) {
  stop(
    "Missing directories: ",
    paste(missing_directories, collapse = ", "),
    ". Run 01_create_project_directories.R first."
  )
}

# Move files one at a time and report what happened.
move_files <- function(files, destination_directory) {
  moved <- 0L

  if (length(files) == 0) {
    message("No files found for: ", destination_directory)
    return(moved)
  }

  for (source_file in files) {
    destination_file <- file.path(
      destination_directory,
      basename(source_file)
    )

    if (file.exists(destination_file)) {
      message("Skipped; destination already exists: ", destination_file)
    } else if (file.rename(source_file, destination_file)) {
      message("Moved: ", source_file, " -> ", destination_file)
      moved <- moved + 1L
    } else {
      warning("Could not move: ", source_file)
    }
  }

  moved
}

# Only top-level files are selected; files already inside folders are untouched.
r_files <- list.files(
  path = ".",
  pattern = "\\.[Rr]$",
  full.names = TRUE,
  recursive = FALSE
)

# Windows cannot move this file while it is running, so copy it instead.
organizer_name <- "02_organize_project_files.R"
organizer_source <- r_files[basename(r_files) == organizer_name]
organizer_destination <- file.path("scripts", organizer_name)
organizer_copied <- 0L

if (length(organizer_source) == 1 && !file.exists(organizer_destination)) {
  if (file.copy(organizer_source, organizer_destination)) {
    message("Copied organizer: ", organizer_destination)
    organizer_copied <- 1L
  } else {
    warning("Could not copy the organizer into the scripts directory.")
  }
}

r_files <- r_files[basename(r_files) != organizer_name]

csv_files <- list.files(
  path = ".",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

plot_files <- list.files(
  path = ".",
  pattern = "\\.(png|jpg|jpeg|pdf|svg|tif|tiff|bmp)$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

moved_scripts <- move_files(r_files, "scripts") + organizer_copied
moved_data <- move_files(csv_files, "data")
moved_plots <- move_files(plot_files, "plots")

cat(
  "Organization complete:",
  moved_scripts, "R files organized,",
  moved_data, "CSV files, and",
  moved_plots, "plot files moved.\n"
)
