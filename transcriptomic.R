library(factoextra)

options(max.print=1000000)

reduced_transcriptomic_df <- read.csv(file="C:/Users/65834/Documents/FYP/gene_afterFtest.csv", header=TRUE, sep=",", row.names="samples") 

reduced_transcriptomic_df <- scale(reduced_transcriptomic_df)

metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)

metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190")), ]
reduced_transcriptomic_df <- reduced_transcriptomic_df[!(row.names(reduced_transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts")), ]




principal_comps <- prcomp(reduced_transcriptomic_df, center = FALSE, scale. = FALSE) 

summary(principal_comps)

standardised_pc_transcriptomic_df_133 <- principal_comps$x

standardised_pc_transcriptomic_df_133 <- standardised_pc_transcriptomic_df_133[,1:31] 

metadata <- as.data.frame(metadata)
standardised_pc_transcriptomic_df_133 <- as.data.frame(standardised_pc_transcriptomic_df_133)

standardised_pc_transcriptomic_df_133 <- cbind(metadata$status, standardised_pc_transcriptomic_df_133)
colnames(standardised_pc_transcriptomic_df_133)[1] <- "status"

standardised_pc_transcriptomic_df_133 <- standardised_pc_transcriptomic_df_133[standardised_pc_transcriptomic_df_133$status != "bipolar disorder diagnosis: Control", ]

standardised_pc_transcriptomic_df_133 <- standardised_pc_transcriptomic_df_133[,-1]




fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, iter.max = 30, nstart = 25, method = "wss", k.max=10)

kmeans_result <- kmeans(standardised_pc_transcriptomic_df_133, 3, nstart = 25)

fviz_cluster(kmeans_result, data = standardised_pc_transcriptomic_df_133[, 1:31], choose.vars=c("PC1", "PC2"))

cluster_allocation <- kmeans_result$cluster 
cluster_allocation <- data.frame(cluster_allocation)

