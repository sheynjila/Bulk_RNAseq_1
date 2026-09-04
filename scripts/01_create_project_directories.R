# ==============================================================================
# CREATE THE PROJECT DIRECTORIES
# Open the RStudio project, then run this script from the project root folder.
# ==============================================================================

directories <- c("scripts", "data", "plots")

for (directory in directories) {
  if (dir.exists(directory)) {
    message("Already exists: ", directory)
  } else {
    dir.create(directory)
    message("Created: ", directory)
  }
}

cat("Project directories are ready.\n")
