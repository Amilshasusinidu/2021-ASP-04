# Load the saved data
saved_data <- readRDS("data/leukemia_checked.rds")
expression_data <- saved_data$expression_data
sample_info <- saved_data$sample_info

# Load the differential expression results
results <- read.csv(
  "results/limma_all_probes.csv",
  stringsAsFactors = FALSE
)

# Check expression values, sample order and group labels
stopifnot(
  identical(
    colnames(expression_data),
    rownames(sample_info)
  ),
  all(is.finite(expression_data)),
  !anyNA(sample_info$group),
  all(sample_info$group %in% c("B-cell ALL", "T-cell ALL"))
)

dir.create("figures", showWarnings = FALSE)

# Identify significant probe sets
significant <- !is.na(results$adj.P.Val) &
  results$adj.P.Val < 0.05

# Assign colours by the direction of the difference
point_colours <- rep("#B0B0B0", nrow(results))
point_colours[significant & results$logFC < 0] <- "#2878B5"
point_colours[significant & results$logFC > 0] <- "#E87524"

# Create the volcano plot
png(
  "figures/volcano_plot.png",
  width = 1600,
  height = 1200,
  res = 180
)

plot(
  results$logFC,
  -log10(pmax(results$adj.P.Val, .Machine$double.xmin)),
  col = point_colours,
  pch = 19,
  cex = 0.5,
  xlab = "Expression difference (T-cell minus B-cell)",
  ylab = "-log10(BH-adjusted p-value)",
  main = "Differential Expression Between T-cell and B-cell ALL"
)

abline(h = -log10(0.05), lty = 2)
abline(v = 0, lty = 3)

legend(
  "topright",
  legend = c(
    "Higher in B-cell ALL",
    "Higher in T-cell ALL",
    "Adjusted p >= 0.05"
  ),
  col = c("#2878B5", "#E87524", "#B0B0B0"),
  pch = 19,
  bty = "n",
  cex = 0.8
)

dev.off()

# Rank probe sets by adjusted p-value
ranked_results <- results[
  order(results$adj.P.Val, results$P.Value),
]

# Select the top 20 significant probe sets
top_probes <- head(
  ranked_results$probe_id[
    !is.na(ranked_results$adj.P.Val) &
      ranked_results$adj.P.Val < 0.05
  ],
  20
)

stopifnot(
  length(top_probes) >= 2,
  all(top_probes %in% rownames(expression_data))
)

# Standardize each selected probe across samples
heatmap_data <- expression_data[top_probes, , drop = FALSE]

stopifnot(all(apply(heatmap_data, 1, sd) > 0))

heatmap_scaled <- t(scale(t(heatmap_data)))

# Define the sample group colours
group_colours <- c(
  "B-cell ALL" = "#2878B5",
  "T-cell ALL" = "#E87524"
)

sample_colours <- unname(
  group_colours[as.character(sample_info$group)]
)

# Create the heatmap
png(
  "figures/top20_heatmap.png",
  width = 2000,
  height = 1600,
  res = 180
)

# Reserve space for the title and colour explanation
par(oma = c(3, 0, 4, 0))

heatmap(
  heatmap_scaled,
  scale = "none",
  col = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(101),
  ColSideColors = sample_colours,
  labCol = rep("", ncol(heatmap_scaled)),
  margins = c(5, 10),
  main = ""
)

# Add the title outside the heatmap panels
mtext(
  "Top 20 Differentially Expressed Probe Sets",
  side = 3,
  outer = TRUE,
  line = 1,
  cex = 1.2,
  font = 2
)

# Explain the standardized expression colours
mtext(
  "Row z-scores: blue = lower, white = average, red = higher expression",
  side = 1,
  outer = TRUE,
  line = 1,
  cex = 0.8
)

# Explain the sample annotation bar
legend(
  "topright",
  legend = names(group_colours),
  fill = group_colours,
  bty = "n",
  cex = 0.7
)

dev.off()

# Display the result summary
cat("\nSignificant probe sets:\n")
print(sum(significant))

cat("\nNew figures:\n")
print(c(
  "figures/volcano_plot.png",
  "figures/top20_heatmap.png"
))