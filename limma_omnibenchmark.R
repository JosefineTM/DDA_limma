#!/usr/bin/env Rscript

# Limma wrapper for OmniBenchmark
# Adapts the existing limma_method.r to OmniBenchmark's interface

library(limma)
library(dplyr)
library(optparse)

# Parse OmniBenchmark arguments
option_list <- list(
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--name", type="character", help="Dataset name"),
  make_option("--data.matrix", type="character", help="Input data file"),
  make_option("--data.true_labels", type="character", help="True labels")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Debug: print all arguments
cat("Running limma analysis\n")
cat("All arguments:\n")
print(args)
cat("\n")

# Access data.matrix with the exact name (including dot)
data_matrix_file <- args$`data.matrix`
cat("Input file:", data_matrix_file, "\n")
cat("Output dir:", args$output_dir, "\n")

# Create output directory
dir.create(args$output_dir, showWarnings=FALSE, recursive=TRUE)

# Read data
data <- read.csv(data_matrix_file, row.names=1)

# Remove label column if present
if ("is_differentially_expressed" %in% colnames(data)) {
  data <- data[, !colnames(data) %in% "is_differentially_expressed"]
}

cat("Data dimensions:", dim(data), "\n")

# Define groups from column names (A* vs B*)
groups <- factor(ifelse(grepl("^A", colnames(data)), "A", "B"))
cat("Groups:", as.character(groups), "\n")

# Create design matrix (without intercept for direct comparison)
design <- model.matrix(~0 + groups)
colnames(design) <- c("groupsA", "groupsB")

# Fit linear model
fit <- lmFit(data, design)

# Create contrast matrix (B vs A)
contrasts_matrix <- makeContrasts(groupsB - groupsA, levels = design)
fit2 <- contrasts.fit(fit, contrasts_matrix)

# Empirical Bayes moderation
fit3 <- eBayes(fit2, robust=TRUE)

# Extract results
results <- topTable(fit3, coef=1, adjust.method="BH", number=Inf, sort.by="none")

# Format output
output_df <- data.frame(
  Name = rownames(data),
  Effect.Size = results$logFC,
  P.Value = results$P.Value,
  adj.P.Val = results$adj.P.Val,
  Significant = results$adj.P.Val < 0.05
)

# Save results
output_file <- file.path(args$output_dir, paste0(args$name, "_limma_results.csv"))
# output_file <- args$output_dir
write.csv(output_df, output_file, row.names=FALSE)


cat("Results saved to:", output_file, "\n")
cat("Total features:", nrow(output_df), "\n")
cat("Significant (FDR < 0.05):", sum(output_df$Significant), "\n")
