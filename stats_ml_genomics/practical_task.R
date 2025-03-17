library(ggplot2)
library(stats)
library(cluster)
library(tidyverse)


# Linear model fitting ----

# Toy dataset with rows representing genes  
# and columns representing histone modification scores (H3K4me3, H3K27me3) and gene expression values (measured_log2)
d <- read.table("./datasets/HistoneModeVSgeneExp.txt", header = TRUE)

# Fit linear models
mod1=lm(d$measured_log2~d$H3k4me3)
summary(mod1)

coefficients1 <- coef(mod1)
print(coefficients1)


mod2=lm(d$measured_log2~d$H3k27me3)
summary(mod2)

coefficients2 <- coef(mod2)
print(coefficients2)


# Create the scatter plot with fitted regression line
ggplot(d, aes(x = H3k4me3, y = measured_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "H3k4me3", 
       y = "Gene Expression (log2)",
       title = "Gene Expression vs H3k4me3 Levels with Fitted Regression Line") +
  theme_minimal()

ggplot(d, aes(x = H3k27me3, y = measured_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "H3K27me3", 
       y = "Gene Expression (log2)",
       title = "Gene Expression vs H3K27me3 Levels with Fitted Regression Line") +
  theme_minimal()



# k-means clustering ----

data <- read.table("./datasets/leukemiaExp.txt", header = TRUE)

# 1. Read the file
gene_names <- rownames(data)  # Save gene names
data <- as.matrix(data)  # Convert to matrix

# 2. Calculate variance for each gene
# The apply function returnes a vector containing the variance for each row.
gene_variances <- apply(data, 1, var)

# 3. Select top 1000 genes with highest variance
top_1000_genes <- order(gene_variances, decreasing = TRUE)[1:1000]
data_filtered <- data[top_1000_genes, ]

# Transpose the data for k-means clustering
data_filtered_t <- t(data_filtered)

# take logs: reduce the effect of highly expressed genes and make the data more "normally distributed"
data_filtered_t_log <- log2(data_filtered_t+1)

# scale: ensure that all genes contribute equally to the clustering (preventing highly expressed genes from dominating)
data_filtered_t_log_sc <- scale(data_filtered_t_log)

# 4 & 5. Perform k-means clustering and calculate Silhouette scores
max_k <- 10  # Maximum number of clusters to try
silhouette_scores <- numeric(max_k - 1)

for (k in 2:max_k) {
  km <- kmeans(data_filtered_t_log_sc, centers = k, nstart = 25)
  silhouette_scores[k-1] <- mean(silhouette(km$cluster, dist(data_filtered_t_log_sc))[, 3])
}

# Plot Silhouette scores
plot(2:max_k, silhouette_scores, type = "b", xlab = "Number of clusters (k)", 
     ylab = "Average Silhouette score", main = "Silhouette Analysis")

# 6. Determine optimum number of clusters
optimum_k <- which.max(silhouette_scores) + 1
cat("The optimum number of clusters is:", optimum_k, "\n")

# Perform final clustering with optimum k
final_clustering <- kmeans(data_filtered_t_log_sc, centers = optimum_k, nstart = 25)

# Print cluster sizes
print(table(final_clustering$cluster))

# 6. Create PCA plot
pca_result <- prcomp(data_filtered_t_log_sc)
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$cluster <- as.factor(final_clustering$cluster)

ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Leukemia Patients",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"))


# Split the rownames (patient IDs) into a new column with the patient subtype info
pca_data$subtype <- sapply(rownames(pca_data), function(x) strsplit(x, "\\.")[[1]][2])
unique(pca_data$subtype)




# Random Foreests ----

library(randomForest)

# Read the data
data <- read.table("./datasets/CpGMeth2Age.txt", header = TRUE)

# Separate features (X) and target variable (y)
X <- as.matrix(data[, -1])  # All columns except the first one
y <- data[[1]]  # First column (age)

# Split the data into training and testing sets
set.seed(42)  # for reproducibility
train_indices <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# Build the random forest model
rf_model <- randomForest(x = X_train, y = y_train, ntree = 500)
summary(rf_model)
print(rf_model)

# Make predictions on the test set
predictions <- predict(rf_model, X_test)

# Calculate R-squared
SSR <- sum((predictions - y_test)^2)
SST <- sum((y_test - mean(y_test))^2)
R_squared <- 1 - SSR/SST

# Print the R-squared value
cat("R-squared:", R_squared, "\n")

# Plot predicted vs actual values
plot(y_test, predictions, main = "Predicted vs Actual Age",
     xlab = "Actual Age", ylab = "Predicted Age")
abline(0, 1, col = "red")
