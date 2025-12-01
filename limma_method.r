options(repos = c(CRAN = "https://cloud.r-project.org"))
# TO REMOVE ONCE WE SETUP THE ENVIRONMENT
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
install.packages("dplyr")
install.packages("tibble")
install.packages("optparse")

library(BiocManager)
library(limma)
library(dplyr)
library(tibble)
library(optparse)

option_list <- list(
  make_option(c("-d", "--dataset")),
  make_option(c("-l", "--label_samples")),
  make_option(c("-m", "--metadata")),
  make_option(c("-o", "--output"))
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# load the dataframe containing the group labels for each sample id
label_samples_df <- read.csv(args$label_samples,row.names = 1)
label_samples_df <- label_samples_df[ , !grepl("^X", names(label_samples_df)) ]
colnames(label_samples_df) <- c("Sample_ID", "Group_Label")
# load the abundance matrix
abundance_matrix <- read.csv(args$dataset,row.names = 1)
colnames(abundance_matrix) <- sub("^X", "", colnames(abundance_matrix))
# load the dataframe containing the gene names and protein IDs
meta_data_df <- read.csv(args$metadata, row.names = 1)
protein_names <- meta_data_df$Name

# --- create group containing the explanotory variable for the DEA analysis
# -- the variables are labels for tumor and normal samples
labels_vec <- unlist(label_samples_df[,2])
groups <- factor(labels_vec)
# create the design matrix. In limma, the design matrix defines the design of the model and contains the information about the parameters to take into account for the linear model to be fit on the data later on
# here we specify a design where parameters to consider are related to the two groups (T vs N) of samples-without the intercept (meaning that we consider the means of the two groups independetly to be compared for DEA)
# we are saying: "we want to see if the gene is differentially expressed in one of the two groups. To do so, we want to compare for a each gene the mean across normal samples with the mean across tumor samples and calculate the difference"
design <- model.matrix(~0 + groups)
# fit the linear model to the expression matrix data of each protein/gene
fit <- lmFit(abundance_matrix, design) 

# Now that the model has estimated the groups means for each gene/protein in the dataset- we can proceed with the contrast matrix
# the contrast matrix is used when no intercept is defined in the design matrix meaning that the linear model is not calculating the differences across the two groups (no reference given)- but just estimating the two means
# to actually get information on these estimation, we need to define a rule that the model can use to compare the two groups.
# For example, in this simple case we are just interested in computing the difference in mean expression between tumor and normal samples
contrasts_matrix <- makeContrasts(groupsT - groupsN, levels = design)
# fit the contrast matrix
fit2 <- contrasts.fit(fit, contrasts_matrix)

# --- Empirical Bayes moderated t-test
# -- used to improve the accuracy of the variance estimates for each feature (gene/protein)
# -- this is usually done when the number of sample is low and the estimated variability for a feature can be the result of random fluctuations or instability due to small sample size
# -- the method compute a prior distribution of the variances of all the features and then shrink each individual variance towards the estimated prior distribtion
fit3 <- eBayes(fit2, robust=TRUE)

# get limma results in the form of a table and reassign to each row the corresponding gene name
limma_results <- topTable(fit3, coef = 1, adjust.method="BH", number=Inf, sort.by = "none") 
limma_results <- limma_results %>%
  mutate(Name = protein_names) %>%   # map names based on original row order
  select(Name, everything()) 

# manipulate the limma table to create the desired output format
limma_results_output <- limma_results[,c("Name", "logFC", "P.Value")]
names(limma_results_output)[names(limma_results_output)=="logFC"] <- "Effect.Size"

# save the output file
write.csv(limma_results_output, args$output, row.names = FALSE)