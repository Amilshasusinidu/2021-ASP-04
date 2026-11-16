# Run from the project root containing README.md, code and data
stopifnot(file.exists("data/leukemia_checked.rds"))

# Install limma if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

library(limma)

# Load the saved expression data and sample information
saved_data <- readRDS("data/leukemia_checked.rds")
expression_data <- saved_data$expression_data
sample_info <- saved_data$sample_info

# Check data quality and sample alignment
stopifnot(
  all(is.finite(expression_data)),
  identical(
    colnames(expression_data),
    rownames(sample_info)
  )
)

group <- factor(
  sample_info$group,
  levels = c("B-cell ALL", "T-cell ALL")
)

stopifnot(!anyNA(group), all(table(group) > 0))

# Create the model for the two groups
design <- model.matrix(~ 0 + group)
colnames(design) <- c("B_cell", "T_cell")

# Compare T-cell ALL against B-cell ALL
contrast <- makeContrasts(
  T_vs_B = T_cell - B_cell,
  levels = design
)

# Fit the linear model and moderate the statistics
fit <- lmFit(expression_data, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Adjust p-values using Benjamini-Hochberg
results <- topTable(
  fit,
  coef = "T_vs_B",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Preserve the probe identifiers
results <- data.frame(
  probe_id = rownames(results),
  results,
  row.names = NULL
)

# Select features with an adjusted p-value below 0.05
significant_results <- results[
  !is.na(results$adj.P.Val) & results$adj.P.Val < 0.05,
]

# Save the results
dir.create("results", showWarnings = FALSE)

write.csv(
  results,
  "results/limma_all_probes.csv",
  row.names = FALSE
)

write.csv(
  significant_results,
  "results/limma_significant_probes.csv",
  row.names = FALSE
)

# Record package versions
writeLines(
  capture.output(sessionInfo()),
  "results/session_info.txt"
)

# Print a summary
cat("\nSample counts:\n")
print(table(group))

cat("\nTotal probe sets tested:\n")
print(nrow(results))

cat("\nSignificant probe sets at adjusted p < 0.05:\n")
print(nrow(significant_results))

cat("\nTop five results:\n")
print(head(results, 5))