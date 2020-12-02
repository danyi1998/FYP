library("DESeq2")

options(max.print=1000000)

transcriptomic_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")
transcriptomic_df <- t(transcriptomic_df)
vst_input <- transcriptomic_df 
vst_input <- t(vst_input)
transcriptomic_df <- as.data.frame(transcriptomic_df)

# remove faulty samples as in GSE124326_BD1 page 2
vst_result <- varianceStabilizingTransformation(vst_input, blind = TRUE)
vst_result <- t(vst_result)
variances <- apply(X=vst_result, MARGIN=2, FUN=var)
sorted_variances <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:500] 
genes_with_top500_variances <- vst_result[, sorted_variances]
pca_result <- prcomp(genes_with_top500_variances, center = TRUE, scale. = TRUE)
pca_result <- pca_result$x
pca_result <- as.data.frame(pca_result)
pc1 <- pca_result[,1]
meta_data <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
meta_data <- t(meta_data)
meta_data <- as.data.frame(meta_data)
pc1_and_metadata <- cbind(pc1, meta_data)
pc1_and_metadata_controls_only <- pc1_and_metadata[pc1_and_metadata$status == "bipolar disorder diagnosis: Control", ]
pc1_and_metadata_controls_only <- pc1_and_metadata_controls_only[!(row.names(pc1_and_metadata_controls_only) %in% c("X12R1998", "X13R1190")), ]

transcriptomic_df <- transcriptomic_df[!(row.names(transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts", "X108773A.counts")), ]
metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)
metadata <- as.data.frame(metadata) 
metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190", "X108773A")), ]



transcriptomic_df <- transcriptomic_df[vapply(transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 
variances <- apply(X=transcriptomic_df, MARGIN=2, FUN=var)
variances <- as.data.frame(variances)
dummy <- vector(length=55659)
dummy <- as.data.frame(dummy)
variances <- cbind(variances, dummy)
variances <- variances[variances$variances > 5, ]

transcriptomic_df <- scale(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

transcriptomic_df <- transcriptomic_df[,row.names(variances)]




# for 3k+ genes
reduced_transcriptomic_df <- read.csv(file="C:/Users/65834/Documents/FYP/gene_afterFtest.csv", header=TRUE, sep=",", row.names="samples") 
reduced_transcriptomic_df <- reduced_transcriptomic_df[!(row.names(reduced_transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts", "X108773A.counts")), ]

metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)
metadata <- as.data.frame(metadata) 
metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190", "X108773A")), ]

reduced_transcriptomic_df <- reduced_transcriptomic_df[vapply(reduced_transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 

reduced_transcriptomic_df <- scale(reduced_transcriptomic_df)
reduced_transcriptomic_df <- as.data.frame(reduced_transcriptomic_df)





principal_comps <- prcomp(reduced_transcriptomic_df, center = FALSE, scale. = FALSE) 

summary(principal_comps)

gene_loadings <- principal_comps$rotation
gene_loadings <- gene_loadings[, 1:133]

standardised_pc_transcriptomic_df_133 <- principal_comps$x
standardised_pc_transcriptomic_df_133 <- as.data.frame(pc_transcriptomic_df)

standardised_pc_transcriptomic_df_133 <- pc_transcriptomic_df[,1:113] 

#standardised_pc_transcriptomic_df_133 <- scale(pc_transcriptomic_df_133)
#standardised_pc_transcriptomic_df_133 <- as.data.frame(standardised_pc_transcriptomic_df_133) 

standardised_pc_transcriptomic_df_133 <- cbind(metadata$status, standardised_pc_transcriptomic_df_133)
colnames(standardised_pc_transcriptomic_df_133)[1] <- "status"

standardised_pc_transcriptomic_df_133 <- standardised_pc_transcriptomic_df_133[standardised_pc_transcriptomic_df_133$status != "bipolar disorder diagnosis: Control", ]

standardised_pc_transcriptomic_df_133 <- standardised_pc_transcriptomic_df_133[,-1]


# hierarchical clustering
fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, method = "wss", k.max=20)

fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, method = "silhouette", k.max=20)

distance <- dist(standardised_pc_transcriptomic_df_133)

hc <- hclust(distance, method = "complete")

plot(hc, labels = FALSE)

rect.hclust(hc, k=6)

cluster_groups <- cutree(hc, k=2)

hclust_allocation <- cbind(standardised_pc_transcriptomic_df_133, cluster_groups)

cluster_groups <- as.data.frame(cluster_groups)




# k means 
library(factoextra)
library(cluster)

set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, iter.max = 30, nstart = 25, method = "wss", k.max=20)
set.seed(123)
fviz_nbclust(standardised_pc_transcriptomic_df_133, kmeans, iter.max = 30, nstart = 25, method = "silhouette", k.max=20)

set.seed(123)
kmeans_result <- kmeans(standardised_pc_transcriptomic_df_133, 3, nstart = 25)

fviz_cluster(kmeans_result, data = standardised_pc_transcriptomic_df_133[, 1:113], choose.vars=c("PC1", "PC2"))

cluster_allocation <- kmeans_result$cluster 
cluster_allocation <- data.frame(cluster_allocation)

standardised_pc_transcriptomic_df_133_extracolumns <- cbind(standardised_pc_transcriptomic_df_133, cluster_allocation)



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

cluster_age <- plot(standardised_pc_transcriptomic_df_133_extracolumns$cluster_allocation, standardised_pc_transcriptomic_df_133_extracolumns$age, type="p", xlab="cluster", ylab="age", xlim=c(0,12))

cluster_sex <- aggregate(standardised_pc_transcriptomic_df_133_extracolumns, by=list(standardised_pc_transcriptomic_df_133_extracolumns$cluster_allocation, standardised_pc_transcriptomic_df_133_extracolumns$sex), FUN=length)
cluster_sex <- cluster_sex[, 1:3]
colnames(cluster_sex) <- c("cluster","sex", "count")
cluster_sex_plot <- plot(cluster_sex$cluster, cluster_sex$count, type="p", xlab="cluster", ylab="count", main="Cluster-Sex" , col=ifelse(cluster_sex$sex=="F",'red','blue'), xlim=c(0,12))

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
