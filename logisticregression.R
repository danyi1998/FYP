library(caret)

options(max.print=1000000)

# classification with pca

# pre-processing
reduced_transcriptomic_df <- read.csv(file="C:/Users/65834/Documents/FYP/gene_afterFtest.csv", header=TRUE, sep=",", row.names="samples") 

reduced_transcriptomic_df <- scale(reduced_transcriptomic_df)

metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)

metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190")), ]
reduced_transcriptomic_df <- reduced_transcriptomic_df[!(row.names(reduced_transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts")), ]


# pca
principal_comps <- prcomp(reduced_transcriptomic_df, center = FALSE, scale. = FALSE) 

standardised_pc_transcriptomic_df_133 <- principal_comps$x
standardised_pc_transcriptomic_df_133 <- standardised_pc_transcriptomic_df_133[,1:31] 

metadata <- as.data.frame(metadata)
standardised_pc_transcriptomic_df_133 <- as.data.frame(standardised_pc_transcriptomic_df_133)

standardised_pc_transcriptomic_df_133 <- cbind(metadata$status, standardised_pc_transcriptomic_df_133)
colnames(standardised_pc_transcriptomic_df_133)[1] <- "status"
standardised_pc_transcriptomic_df_133$status <- gsub("bipolar disorder diagnosis: ", "", standardised_pc_transcriptomic_df_133$status)
standardised_pc_transcriptomic_df_133$status <- gsub("BP1", "Patient", standardised_pc_transcriptomic_df_133$status)
standardised_pc_transcriptomic_df_133$status <- gsub("BP2", "Patient", standardised_pc_transcriptomic_df_133$status)


# classification
standardised_pc_transcriptomic_df_133$status <- as.factor(standardised_pc_transcriptomic_df_133$status) 

cross_validation <- trainControl(method = "cv", number = 5)

model <- train(status~., data=standardised_pc_transcriptomic_df_133, method="glm", family="binomial", maxit=40, trControl=cross_validation)

print(model)









# classification without pca

# pre-processing
reduced_transcriptomic_df <- read.csv(file="C:/Users/65834/Documents/FYP/gene_afterFtest.csv", header=TRUE, sep=",", row.names="samples") 

reduced_transcriptomic_df <- scale(reduced_transcriptomic_df)

metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)

metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190")), ]
reduced_transcriptomic_df <- reduced_transcriptomic_df[!(row.names(reduced_transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts")), ]


# preparation 
metadata <- as.data.frame(metadata)
reduced_transcriptomic_df <- as.data.frame(reduced_transcriptomic_df)

reduced_transcriptomic_df <- cbind(metadata$status, reduced_transcriptomic_df)
colnames(reduced_transcriptomic_df)[1] <- "status"
reduced_transcriptomic_df$status <- gsub("bipolar disorder diagnosis: ", "", reduced_transcriptomic_df$status)
reduced_transcriptomic_df$status <- gsub("BP1", "Patient", reduced_transcriptomic_df$status)
reduced_transcriptomic_df$status <- gsub("BP2", "Patient", reduced_transcriptomic_df$status)


# classification
reduced_transcriptomic_df$status <- as.factor(reduced_transcriptomic_df$status) 

cross_validation <- trainControl(method = "cv", number = 5)

model <- train(status~., data=reduced_transcriptomic_df, method="glm", family="binomial", maxit=40, trControl=cross_validation)

print(model)



