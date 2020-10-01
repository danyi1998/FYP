transcriptomic_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")
transcriptomic_df <- t(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

transcriptomic_df <- transcriptomic_df[vapply(transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 

principal_comps <- prcomp(transcriptomic_df, center = TRUE, scale. = TRUE) 

summary(principal_comps)

pc_transcriptomic_df <- principal_comps$x
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

pc_transcriptomic_df_2 <- pc_transcriptomic_df[,1:2] 

standardised_pc_transcriptomic_df_2 <- scale(pc_transcriptomic_df_2)
standardised_pc_transcriptomic_df_2 <- as.data.frame(standardised_pc_transcriptomic_df_2)

library(factoextra)
library(cluster)

set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_2, kmeans, iter.max = 30, nstart = 25, method = "wss")
set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_2, kmeans, iter.max = 30, nstart = 25, method = "silhouette")

set.seed(123)
kmeans_result <- kmeans(standardised_pc_transcriptomic_df_2, 5, nstart = 25)

fviz_cluster(kmeans_result, data = standardised_pc_transcriptomic_df_2[, 1:2])

cluster_allocation <- kmeans_result$cluster 
cluster_allocation <- data.frame(cluster_allocation)

standardised_pc_transcriptomic_df_2_extracolumns <- cbind(standardised_pc_transcriptomic_df_2, cluster_allocation)



metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)
metadata <- as.data.frame(metadata)

metadata$status <- gsub("bipolar disorder diagnosis: ", "", metadata$status)
metadata$age <- gsub("age: ", "", metadata$age)
metadata$sex <- gsub("Sex: ", "", metadata$sex)
metadata$lithium <- gsub("lithium use", "", metadata$lithium)
metadata$lithium <- gsub("non-user=0, user = 1", "", metadata$lithium)
metadata$lithium <- gsub("[ ():]", "", metadata$lithium)
metadata$tobacco <- gsub("tobacco use: ", "", metadata$tobacco)

metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190")), ]
standardised_pc_transcriptomic_df_2_extracolumns <- standardised_pc_transcriptomic_df_2_extracolumns[!(row.names(standardised_pc_transcriptomic_df_2_extracolumns) %in% c("X12R1998.counts", "X13R1190.counts")), ]

metadata$age <- as.numeric(metadata$age)

standardised_pc_transcriptomic_df_2_extracolumns <- cbind(standardised_pc_transcriptomic_df_2_extracolumns, metadata)

standardised_pc_transcriptomic_df_2_extracolumns <- standardised_pc_transcriptomic_df_2_extracolumns[standardised_pc_transcriptomic_df_2_extracolumns$status != "Control", ]
