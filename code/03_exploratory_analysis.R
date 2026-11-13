# Run 01_import_data.R and 02_check_groups.R first

# Check the data and sample order
stopifnot(all(is.finite(expression_data)))
stopifnot(
  identical(
    colnames(expression_data),
    rownames(sample_info)
  )
)
stopifnot(!anyNA(sample_info$group))

# Create a folder for the figures
dir.create("figures", showWarnings = FALSE)

# Set colours for the two groups
group_colours <- c(
  "B-cell ALL" = "#2878B5",
  "T-cell ALL" = "#E87524"
)

sample_colours <- unname(
  group_colours[as.character(sample_info$group)]
)

# Save a boxplot of all samples
png(
  "figures/expression_boxplot.png",
  width = 1800,
  height = 1000,
  res = 150
)

boxplot(
  expression_data,
  col = sample_colours,
  border = sample_colours,
  outline = FALSE,
  names = FALSE,
  main = "Expression Distributions Across Samples",
  xlab = "Samples in dataset order",
  ylab = "Processed expression value"
)

legend(
  "topright",
  legend = names(group_colours),
  fill = group_colours,
  bty = "n"
)

dev.off()

# Keep expression features that vary across samples
feature_variance <- apply(expression_data, 1, var)
pca_data <- expression_data[feature_variance > 0, , drop = FALSE]

# PCA requires samples as rows and features as columns
pca_result <- prcomp(
  t(pca_data),
  center = TRUE,
  scale. = FALSE
)

variance_percent <- (
  pca_result$sdev^2 / sum(pca_result$sdev^2)
) * 100

# Save the PCA plot
png(
  "figures/pca_plot.png",
  width = 1400,
  height = 1000,
  res = 150
)

plot(
  pca_result$x[, 1],
  pca_result$x[, 2],
  col = sample_colours,
  pch = 19,
  xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
  ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
  main = "PCA of Childhood Leukemia Samples"
)

legend(
  "topright",
  legend = names(group_colours),
  col = group_colours,
  pch = 19,
  bty = "n"
)

dev.off()

# Print a short summary
print(table(sample_info$group))
print(round(variance_percent[1:2], 2))

cat("\nFigures saved in:\n")
print(normalizePath("figures"))