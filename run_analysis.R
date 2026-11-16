# Set the working directory to this repository root before running.
stopifnot(file.exists("README.md"), dir.exists("code"))
source("code/01_import_data.R")
source("code/02_check_groups.R")
source("code/03_exploratory_analysis.R")
source("code/04_differential_expression.R")
source("code/05_results_visualisation.R")
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
