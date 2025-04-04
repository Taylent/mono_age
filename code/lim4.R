library(limma)
library(dplyr)
library(lme4)
library(car)
install.packages("preprocessCore")
library(preprocessCore)

gene_data <- read.csv("D:/Desktop/data.csv")

# Suppose you have loaded the data into a dataframe named "gene_data"

# Normalization and log2 transformation
normalized_data <- log2(gene_data[, -c(1:4)] + 1)

age <- gene_data$age
gene_data <- na.omit(gene_data)
normalized_data <- normalized_data[1:nrow(gene_data),]
# Fit a mixed-effects model
sum(!complete.cases(gene_data))  # Calculate the number of missing values in all rows of gene_data
response_vector <- rowMeans(normalized_data)

# Then create a dataframe containing the new response variable
new_gene_expr_data <- gene_data  # Copy the original dataframe
new_gene_expr_data$normalized_response <- response_vector  # Add a new response variable column
# Calculate the correlation coefficient between each gene and age
correlation <- apply(normalized_data[,-1], 2, function(gene) cor(age, gene))

# Perform significance testing and extract significantly correlated genes
significance_threshold <- 0.05
significant_genes <- character(length(correlation))  # Create an empty character vector to store the names of significant genes
for (i in 1:length(correlation)) {
  if (sd(normalized_data[, i]) == 0) {
    next  # If the standard deviation is zero, skip this gene
  }
  p_value <- cor.test(age, normalized_data[, i])$p.value
  if (p_value <= significance_threshold) {
    significant_genes[i] <- names(normalized_data)[i]
  }
}
significant_genes <- significant_genes[!is.na(significant_genes)]  # Remove null values

significant_genes
model_list <- list()

# Loop through each gene
for (gene in colnames(normalized_data)) {
  # Build a linear model
  model <- lm(paste(gene, "~ age + gender+ (1 | group)"), data = new_gene_expr_data)
  
  # Store the model results in the list
  model_list[[gene]] <- summary(model)
}

# View the summary of the first gene's model
summary(model_list[[2]])

# Calculate the magnitude of fixed effects (assuming model results are stored in model_list)
fixed_effects <- lapply(model_list, function(model) coef(model))
fixed_effects
# Create an empty dataframe to store results
results <- lapply(model_list, function(model_summary) {
  parameter_names <- row.names(model_summary$coefficients)
  p_values <- model_summary$coefficients[, "Pr(>|t|)"]
  estimates <- model_summary$coefficients[, "Estimate"]
  
  data.frame(
    model = rep(names(model_summary), length(p_values)),
    parameter = parameter_names,
    p_value = p_values,
    estimate = estimates,
    stringsAsFactors = FALSE
  )
})

# Combine all results into one dataframe
results_combined <- do.call(rbind, results)

# Perform FDR correction on p-values
results_combined$fdr <- p.adjust(results_combined$p_value, method = "fdr")

# Select significant results at the 0.05 significance level
significant_results <- subset(results_combined, fdr < 0.05)

# Print the results
print(significant_results)
