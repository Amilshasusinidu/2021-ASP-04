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

# Download the gene expression dataset
geo_data <- getGEO(
  "GSE13425",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# Display the available platforms
print(names(geo_data))

# Display dataset dimensions and sample information
for (dataset in geo_data) {
  print(dim(exprs(dataset)))
  print(head(pData(dataset)))
}

