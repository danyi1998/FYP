options(max.print=1000000)

transcriptomic_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")
transcriptomic_df <- t(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

transcriptomic_df <- transcriptomic_df[vapply(transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 
transcriptomic_df <- scale(transcriptomic_df)

principal_comps <- prcomp(transcriptomic_df, center = TRUE, scale. = TRUE) 

summary(principal_comps)

gene_loadings <- principal_comps$rotation
gene_loadings <- gene_loadings[, 1:133]

pc_transcriptomic_df <- principal_comps$x
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

pc_transcriptomic_df_133 <- pc_transcriptomic_df[,1:133] 

standardised_pc_transcriptomic_df_133 <- scale(pc_transcriptomic_df_133)
standardised_pc_transcriptomic_df_133 <- as.data.frame(standardised_pc_transcriptomic_df_133)

library(factoextra)
library(cluster)

set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, iter.max = 30, nstart = 25, method = "wss")
set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, iter.max = 30, nstart = 25, method = "silhouette")

set.seed(123)
kmeans_result <- kmeans(standardised_pc_transcriptomic_df_133, 5, nstart = 25)

fviz_cluster(kmeans_result, data = standardised_pc_transcriptomic_df_133[, 1:133])

cluster_allocation <- kmeans_result$cluster 
cluster_allocation <- data.frame(cluster_allocation)

standardised_pc_transcriptomic_df_133_extracolumns <- cbind(standardised_pc_transcriptomic_df_133, cluster_allocation)



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

metadata_retain_faulty_samples <- metadata

metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190")), ]
standardised_pc_transcriptomic_df_133_extracolumns <- standardised_pc_transcriptomic_df_133_extracolumns[!(row.names(standardised_pc_transcriptomic_df_133_extracolumns) %in% c("X12R1998.counts", "X13R1190.counts")), ]

metadata$age <- as.numeric(metadata$age)

standardised_pc_transcriptomic_df_133_extracolumns <- cbind(standardised_pc_transcriptomic_df_133_extracolumns, metadata)

standardised_pc_transcriptomic_df_133_extracolumns <- standardised_pc_transcriptomic_df_133_extracolumns[standardised_pc_transcriptomic_df_133_extracolumns$status != "Control", ]

cluster_age <- plot(standardised_pc_transcriptomic_df_133_extracolumns$cluster_allocation, standardised_pc_transcriptomic_df_133_extracolumns$age, type="p", xlab="cluster", ylab="age")

cluster_sex <- aggregate(standardised_pc_transcriptomic_df_133_extracolumns, by=list(standardised_pc_transcriptomic_df_133_extracolumns$cluster_allocation, standardised_pc_transcriptomic_df_133_extracolumns$sex), FUN=length)
cluster_sex <- cluster_sex[, 1:3]
colnames(cluster_sex) <- c("cluster","sex", "count")
cluster_sex_plot <- plot(cluster_sex$cluster, cluster_sex$count, type="p", xlab="cluster", ylab="count", main="Cluster-Sex" , col=ifelse(cluster_sex$sex=="F",'red','blue'))

cluster_lithium <- aggregate(standardised_pc_transcriptomic_df_133_extracolumns, by=list(standardised_pc_transcriptomic_df_133_extracolumns$cluster_allocation, standardised_pc_transcriptomic_df_133_extracolumns$lithium), FUN=length)
cluster_lithium <- cluster_lithium[, 1:3]
colnames(cluster_lithium) <- c("cluster","lithium", "count")
cluster_lithium_plot <- plot(cluster_lithium$cluster, cluster_lithium$count, type="p", xlab="cluster", ylab="count", main="Cluster-Lithium" , col=ifelse(cluster_lithium$lithium=="1",'red','blue'))

cluster_tobacco <- aggregate(standardised_pc_transcriptomic_df_133_extracolumns, by=list(standardised_pc_transcriptomic_df_133_extracolumns$cluster_allocation, standardised_pc_transcriptomic_df_133_extracolumns$tobacco), FUN=length)
cluster_tobacco <- cluster_tobacco[, 1:3]
colnames(cluster_tobacco) <- c("cluster","tobacco", "count")
cluster_tobacco_plot <- plot(cluster_tobacco$cluster, cluster_tobacco$count, type="p", xlab="cluster", ylab="count", main="Cluster-Tobacco" , col=ifelse(cluster_tobacco$tobacco=="1",'red','blue'))



transcriptomic_genes_with_cluster_allocation <- cbind(transcriptomic_df, cluster_allocation)
sample_status_retain_faulty_samples <- metadata_retain_faulty_samples[, 1]
sample_status_retain_faulty_samples <- as.data.frame(sample_status_retain_faulty_samples)
transcriptomic_genes_with_cluster_allocation_and_status <- cbind(transcriptomic_genes_with_cluster_allocation, sample_status_retain_faulty_samples)
transcriptomic_genes_with_cluster_allocation_and_status_no_controls <- transcriptomic_genes_with_cluster_allocation_and_status[transcriptomic_genes_with_cluster_allocation_and_status$sample_status_retain_faulty_samples != "Control", ]
high_loading_gene <- transcriptomic_genes_with_cluster_allocation_and_status_no_controls[, c("ENSG00000244161.1", "cluster_allocation")]


# investigating correlation between gene and clinical data, put on hold for now, check again when using in the future 
metadata_retain_faulty_samples_no_controls <- metadata_retain_faulty_samples[metadata_retain_faulty_samples$status != "Control", ]
metadata_retain_faulty_samples_no_controls$age <- as.numeric(metadata_retain_faulty_samples_no_controls$age)
metadata_retain_faulty_samples_no_controls$lithium <- as.numeric(metadata_retain_faulty_samples_no_controls$lithium)
corr_results <- cor(x=high_loading_gene$ENSG00000121390.13, y=metadata_retain_faulty_samples_no_controls$lithium, use="all.obs", method="pearson")
