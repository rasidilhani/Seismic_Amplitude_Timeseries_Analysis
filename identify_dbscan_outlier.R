library(easylabel)
library(dplyr)

# Load your PCA + cluster results (from your DBSCAN/HDBSCAN pipeline)
plot_df <- data.frame(
  Dim1 = pca_features[,1],
  Dim2 = pca_features[,2],
  cluster = factor(dbscan_result$cluster),
  point_id = df_feat$unix_timestamp  # or datetime_nz
)

easylabel(
  plot_df,
  x = 'Dim1',
  y = 'Dim2',
  col = 'cluster',
  labs = 'point_id',
  size = 8
)