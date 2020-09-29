transcriptomic_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")
transcriptomic_df <- t(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

transcriptomic_df <- transcriptomic_df[vapply(transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 

principal_comps <- prcomp(transcriptomic_df, center = TRUE, scale. = TRUE) 

pc_transcriptomic_df <- principal_comps$x
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

pc_transcriptomic_df_10 <- pc_transcriptomic_df[,1:10] 

standardised_pc_transcriptomic_df_10 <- scale(pc_transcriptomic_df_10)
standardised_pc_transcriptomic_df_10 <- as.data.frame(standardised_pc_transcriptomic_df_10)

library(factoextra)
library(cluster)

set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_10, kmeans, iter.max = 30, nstart = 25, method = "wss")
set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_10, kmeans, iter.max = 30, nstart = 25, method = "silhouette")

set.seed(123)
kmeans_result <- kmeans(standardised_pc_transcriptomic_df_10, 8, nstart = 25)

fviz_cluster(kmeans_result, data = standardised_pc_transcriptomic_df_10[, 1:2])

cluster_allocation <- kmeans_result$cluster 
cluster_allocation <- data.frame(cluster_allocation)

standardised_pc_transcriptomic_df_10_extracolumns <- cbind(standardised_pc_transcriptomic_df_10, cluster_allocation)
