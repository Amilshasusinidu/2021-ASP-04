# Load the required package
library(GEOquery)

# Extract expression data and sample information
dataset <- geo_data[[1]]
expression_data <- exprs(dataset)
sample_info <- pData(dataset)

# Read the sample labels
labels <- trimws(
  as.character(sample_info$characteristics_ch1)
)

# Create a group for each sample
group <- rep(NA_character_, length(labels))

# Allow one or more spaces between Precursor-B and ALL
group[grepl("^Precursor-B[[:space:]]+ALL", labels)] <- "B-cell ALL"
group[labels == "T-ALL"] <- "T-cell ALL"

# Check that every sample belongs to a recognised group
stopifnot(!anyNA(group))

# Check that sample order matches the expression data
stopifnot(
  identical(
    colnames(expression_data),
    rownames(sample_info)
  )
)

# Add the group to the sample information
sample_info$group <- factor(
  group,
  levels = c("B-cell ALL", "T-cell ALL")
)

# Display sample counts
cat("\nSample counts:\n")
print(table(sample_info$group))

# Display expression features and sample count
cat("\nExpression matrix dimensions:\n")
print(dim(expression_data))

# Display the number of missing expression values
cat("\nMissing expression values:\n")
print(sum(is.na(expression_data)))
