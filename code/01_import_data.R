# Install the required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

# Load GEOquery
library(GEOquery)

# Allow more time for the download
options(timeout = 600)

# Download the matrix directly, avoiding GEO directory-listing requests.
# Reuse a successfully downloaded file on later runs.
dir.create("data", showWarnings = FALSE)
matrix_file <- "data/GSE13425_series_matrix.txt.gz"

if (!file.exists(matrix_file)) {
  matrix_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    "GSE13nnn/GSE13425/matrix/GSE13425_series_matrix.txt.gz"
  )
  temporary_file <- paste0(matrix_file, ".download")
  download.file(matrix_url, temporary_file, mode = "wb", method = "libcurl")
  stopifnot(file.rename(temporary_file, matrix_file))
}

# A local series matrix returns an ExpressionSet. Preserve the list
# representation expected by the remaining scripts.
geo_data <- list(
  GSE13425 = getGEO(filename = matrix_file, getGPL = FALSE)
)

# Display the available platforms
print(names(geo_data))

# Display dataset dimensions and sample information
for (dataset in geo_data) {
  print(dim(exprs(dataset)))
  print(head(pData(dataset)))
}

