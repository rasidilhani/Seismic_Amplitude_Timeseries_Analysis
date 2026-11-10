library(tidyverse)

set.seed(123456789, kind = "Mersenne-Twister")


library(readr)
data <- read_csv("data/weekly_HC_RSAM_WIZ_20140701_to_20191231.csv")
View(data)


feature_cols <- c(
  "Fisher_Info",
  "C_Fisher",
  "Var_H_Shannon",
  "Var_C_Shannon",
  "Var_H_Tsallis",
  "Var_H_Renyi",
  "Disequilibrium"
)

-----------------------

kmeans_data <- data[,c(10,11,16,17,18,19,26)]
dim(data)
# ------------------------------------------------------------
# 2. Scale the features
# ------------------------------------------------------------

X <- scale(kmeans_data[, feature_cols])


# ------------------------------------------------------------
# 3. Run K-means with 3 clusters
# ------------------------------------------------------------

set.seed(123)

km_model <- kmeans(
  X,
  centers = 3,
  nstart = 50
)
print(km_model)

# ------------------------------------------------------------
# 4. Add cluster number to the original data
# ------------------------------------------------------------

kmeans_data$Kmeans_Cluster <- factor(km_model$cluster)


# ------------------------------------------------------------
# 5. Number of observations in each cluster
# ------------------------------------------------------------

table(kmeans_data$Kmeans_Cluster)


# ------------------------------------------------------------
# 6. Compare K-means clusters with the true VAL classes
# ------------------------------------------------------------

table(
  True_Class = kmeans_data$VAL_Class,
  Kmeans_Cluster = kmeans_data$Kmeans_Cluster
)


# ------------------------------------------------------------
# 7. Cluster centres
# ------------------------------------------------------------

print(km_model$centers)