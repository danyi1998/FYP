library(e1071)
library(caret)
library(kernlab)

options(max.print=1000000)

transcriptomic_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")
transcriptomic_df <- t(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

transcriptomic_df <- transcriptomic_df[!(row.names(transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts", "X108773A.counts")), ]
metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)
metadata <- as.data.frame(metadata) 
metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190", "X108773A")), ]

transcriptomic_df <- transcriptomic_df[vapply(transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 

variances <- apply(X=transcriptomic_df, MARGIN=2, FUN=var)
sorted_variances <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:10000] 
bottom_variance_genes <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[10001:16500] 

transcriptomic_df <- scale(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

genes_with_highest_variances <- transcriptomic_df[, sorted_variances]
genes_with_lowest_variances <- transcriptomic_df[, bottom_variance_genes]
genes_with_selected_variances <- cbind(genes_with_highest_variances, genes_with_lowest_variances)



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





# naive bayes before pca 
sample_status_and_transcriptomic_df <- cbind(metadata$status, reduced_transcriptomic_df)
colnames(sample_status_and_transcriptomic_df)[1] <- "status"
sample_status_and_transcriptomic_df$status <- gsub("bipolar disorder diagnosis: ", "", sample_status_and_transcriptomic_df$status)
sample_status_and_transcriptomic_df$status <- gsub("BP1", "Patient", sample_status_and_transcriptomic_df$status)
sample_status_and_transcriptomic_df$status <- gsub("BP2", "Patient", sample_status_and_transcriptomic_df$status)


s <- sample(nrow(sample_status_and_transcriptomic_df), nrow(sample_status_and_transcriptomic_df)*0.8)
training_set <- sample_status_and_transcriptomic_df[s,]
testing_set <- sample_status_and_transcriptomic_df[-s,]

nb <- naiveBayes(as.factor(status)~., data=training_set)

pred <- predict(nb, testing_set[,-1])

confusionMatrix(table(pred, testing_set$status))





# svm before pca
sample_status_and_transcriptomic_df <- cbind(metadata$status, reduced_transcriptomic_df)
colnames(sample_status_and_transcriptomic_df)[1] <- "status"
sample_status_and_transcriptomic_df$status <- gsub("bipolar disorder diagnosis: ", "", sample_status_and_transcriptomic_df$status)
sample_status_and_transcriptomic_df$status <- gsub("BP1", "Patient", sample_status_and_transcriptomic_df$status)
sample_status_and_transcriptomic_df$status <- gsub("BP2", "Patient", sample_status_and_transcriptomic_df$status)

cross_validation <- trainControl(method = "cv", number = 10)

model <- train(as.factor(status)~., data=sample_status_and_transcriptomic_df, method="svmLinear", trControl=cross_validation)

print(model)





# pca
principal_comps <- prcomp(reduced_transcriptomic_df, center = FALSE, scale. = FALSE) 

pc_transcriptomic_df <- principal_comps$x
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

pc_transcriptomic_df_133 <- pc_transcriptomic_df[,1:113] 

#standardised_pc_transcriptomic_df_133 <- scale(pc_transcriptomic_df_133)
#standardised_pc_transcriptomic_df_133 <- as.data.frame(standardised_pc_transcriptomic_df_133)







# naive bayes after pca 
sample_status_and_standardised_pc_transcriptomic_df_133 <- cbind(metadata$status, standardised_pc_transcriptomic_df_133)
colnames(sample_status_and_standardised_pc_transcriptomic_df_133)[1] <- "status"
sample_status_and_standardised_pc_transcriptomic_df_133$status <- gsub("bipolar disorder diagnosis: ", "", sample_status_and_standardised_pc_transcriptomic_df_133$status)
sample_status_and_standardised_pc_transcriptomic_df_133$status <- gsub("BP1", "Patient", sample_status_and_standardised_pc_transcriptomic_df_133$status)
sample_status_and_standardised_pc_transcriptomic_df_133$status <- gsub("BP2", "Patient", sample_status_and_standardised_pc_transcriptomic_df_133$status)


s <- sample(nrow(sample_status_and_standardised_pc_transcriptomic_df_133), nrow(sample_status_and_standardised_pc_transcriptomic_df_133)*0.8)
training_set <- sample_status_and_standardised_pc_transcriptomic_df_133[s,]
testing_set <- sample_status_and_standardised_pc_transcriptomic_df_133[-s,]

nb <- naiveBayes(as.factor(status)~., data=training_set)

pred <- predict(nb, testing_set[,-1])

confusionMatrix(table(pred, testing_set$status))





# svm after pca
sample_status_and_standardised_pc_transcriptomic_df_133 <- cbind(metadata$status, pc_transcriptomic_df_133)
colnames(sample_status_and_standardised_pc_transcriptomic_df_133)[1] <- "status"
sample_status_and_standardised_pc_transcriptomic_df_133$status <- gsub("bipolar disorder diagnosis: ", "", sample_status_and_standardised_pc_transcriptomic_df_133$status)
sample_status_and_standardised_pc_transcriptomic_df_133$status <- gsub("BP1", "Patient", sample_status_and_standardised_pc_transcriptomic_df_133$status)
sample_status_and_standardised_pc_transcriptomic_df_133$status <- gsub("BP2", "Patient", sample_status_and_standardised_pc_transcriptomic_df_133$status)

cross_validation <- trainControl(method = "cv", number = 10)

model <- train(as.factor(status)~., data=sample_status_and_standardised_pc_transcriptomic_df_133, method="svmLinear", trControl=cross_validation)

print(model)


