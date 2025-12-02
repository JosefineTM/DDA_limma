options(repos = c(CRAN = "https://cloud.r-project.org"))
# TO REMOVE ONCE WE SETUP THE ENVIRONMENT
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("limma")
#install.packages("dplyr")
#install.packages("tibble")
#install.packages("optparse")

library(BiocManager)
library(limma)
library(dplyr)
library(tibble)
library(optparse)

option_list <- list(
  make_option(c("-d", "--data.matrix")),
  make_option(c("-l", "--data.true_labels")),
  make_option(c("-o", "--output_dir")), 
  make_option(c("-n", "--name")) 
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

label_samples_df <- read.csv(args$data.true_labels, header= FALSE)
#label_samples_df <- label_samples_df[ , !grepl("^X", names(label_samples_df)) ]
#colnames(label_samples_df) <- c("Sample_ID", "Group_Label")

# load the abundance matrix
abundance_matrix <- read.csv(args$data.matrix, check.names = FALSE, stringsAsFactors = FALSE)
rownames(abundance_matrix) <- make.unique(abundance_matrix[[1]])
abundance_matrix[[1]] <- NULL
#head(abundance_matrix,2)

# create group containing the explanotory variable for the DEA analysis
# -- the variables are labels for tumor and normal samples
labels_vec <- unlist(label_samples_df[,2])
groups <- factor(labels_vec)
# create the design matrix
design <- model.matrix(~0 + groups)
# fit the linear model to the expression matrix data of each protein/gene
fit <- lmFit(abundance_matrix, design) 

# proceed with the contrast matrix
design_cols <- colnames(design)
contrast_exp <- paste0(design_cols[2], " - ", design_cols[1])
contrasts_matrix <- makeContrasts(contrasts = contrast_exp, levels = design)
# fit the contrast matrix
fit2 <- contrasts.fit(fit, contrasts_matrix)
# empirical Bayes moderated t-test
fit3 <- eBayes(fit2, robust=TRUE)

# get limma results in the form of a table and reassign to each row the corresponding gene name
limma_results <- topTable(fit3, coef = 1, adjust.method="BH", number=Inf, sort.by = "none") 

# manipulate the limma table to create the desired output format
limma_results_output <- limma_results[,c("logFC", "P.Value")]
names(limma_results_output)[names(limma_results_output)=="logFC"] <- "effect_size"
names(limma_results_output)[names(limma_results_output)=="P.Value"] <- "P.Value"
# move rownames to a column named "ID"
limma_results_output$ID <- rownames(limma_results)
# remove row names
rownames(limma_results_output) <- NULL
limma_results_output <- limma_results_output[, c("ID", setdiff(names(limma_results_output), "ID"))]
head(limma_results_output)
# save the output file
# create the directory if it doesn't exist
args$output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir,recursive = TRUE, showWarnings = FALSE)
}
#output_filename <- paste0(args$out_file_name,"_results.csv")
#write.csv(limma_results_output, file = file.path(args$output_dir, output_filename))
print(args)
args$data.matrix
args$data.true_labels
args$output_dir
args$name
#output_file
output_filename <- paste0(args[["name"]], "_results.csv")
write.csv(limma_results_output, file=file.path(args[["output_dir"]], output_filename))
print(paste("Results written to:", output_filename))
print(paste("Output directory:", args$output_dir))
print(paste("Output file:", file.path(args$output_dir, output_filename)))