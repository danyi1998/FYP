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






# relationship with clinical data

metadata_patients_only <- metadata[metadata$status != "bipolar disorder diagnosis: Control", ]

metadata_with_cluster <- cbind(cluster_allocation, metadata_patients_only)


metadata_with_cluster$age <- gsub("age: ", "", metadata_with_cluster$age)
metadata_with_cluster$age <- as.numeric(metadata_with_cluster$age)
cluster_age <- plot(metadata_with_cluster$cluster_allocation, metadata_with_cluster$age, type="p", xlab="cluster", ylab="age", xlim=c(0,5))


cluster_sex <- aggregate(metadata_with_cluster, by=list(metadata_with_cluster$cluster_allocation, metadata_with_cluster$sex), FUN=length)
cluster_sex <- cluster_sex[, 1:3]
colnames(cluster_sex) <- c("cluster","sex", "count")
cluster_sex_plot <- plot(cluster_sex$cluster, cluster_sex$count, type="p", xlab="cluster", ylab="count", main="Cluster-Sex" , col=ifelse(cluster_sex$sex=="Sex: F",'red','blue'), xlim=c(0,5))


cluster_lithium <- aggregate(metadata_with_cluster, by=list(metadata_with_cluster$cluster_allocation, metadata_with_cluster$lithium), FUN=length)
cluster_lithium <- cluster_lithium[, 1:3]
colnames(cluster_lithium) <- c("cluster","lithium", "count")
cluster_lithium_plot <- plot(cluster_lithium$cluster, cluster_lithium$count, type="p", xlab="cluster", ylab="count", main="Cluster-Lithium" , col=ifelse(cluster_lithium$lithium=="lithium use (non-user=0, user = 1): 1",'red','blue'), xlim=c(0,5))


cluster_tobacco <- aggregate(metadata_with_cluster, by=list(metadata_with_cluster$cluster_allocation, metadata_with_cluster$tobacco), FUN=length)
cluster_tobacco <- cluster_tobacco[, 1:3]
colnames(cluster_tobacco) <- c("cluster","tobacco", "count")
cluster_tobacco <- cluster_tobacco[cluster_tobacco$tobacco != "tobacco use: NA", ]
cluster_tobacco_plot <- plot(cluster_tobacco$cluster, cluster_tobacco$count, type="p", xlab="cluster", ylab="count", main="Cluster-Tobacco" , col=ifelse(cluster_tobacco$tobacco=="tobacco use: 1",'red','blue'))


