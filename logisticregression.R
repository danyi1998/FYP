library(caret)

options(max.print=1000000)

reduced_transcriptomic_df <- read.csv(file="C:/Users/65834/Documents/FYP/gene_afterFtest.csv", header=TRUE, sep=",", row.names="samples") 
reduced_transcriptomic_df <- reduced_transcriptomic_df[!(row.names(reduced_transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts", "X108773A.counts")), ]

metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)
metadata <- t(metadata)
metadata <- as.data.frame(metadata) 
metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190", "X108773A")), ]
metadata$status <- gsub("bipolar disorder diagnosis: ", "", metadata$status)
metadata$status <- gsub("BP1", "Patients", metadata$status)
metadata$status <- gsub("BP2", "Patients", metadata$status)

reduced_transcriptomic_df <- reduced_transcriptomic_df[vapply(reduced_transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 

reduced_transcriptomic_df <- scale(reduced_transcriptomic_df)
reduced_transcriptomic_df <- as.data.frame(reduced_transcriptomic_df)

reduced_transcriptomic_df <- cbind(metadata$status, reduced_transcriptomic_df)
colnames(reduced_transcriptomic_df)[1] <- "status" 



# start classification here without pca
reduced_transcriptomic_df$status <- as.factor(reduced_transcriptomic_df$status) 
s <- sample(nrow(reduced_transcriptomic_df), nrow(reduced_transcriptomic_df)*0.8)
training_set <- reduced_transcriptomic_df[s,]
testing_set <- reduced_transcriptomic_df[-s,]

model <- glm(status~., training_set, family="binomial", maxit=100)

res <- predict(model, testing_set[,-1], type="response")

table(actual=testing_set$status, predicted=res>0.5)



# with pca
reduced_transcriptomic_df <- reduced_transcriptomic_df[,-1]

principal_comps <- prcomp(reduced_transcriptomic_df, center = FALSE, scale. = FALSE) 

pc_transcriptomic_df <- principal_comps$x
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

pc_transcriptomic_df <- scale(pc_transcriptomic_df)
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

pc_transcriptomic_df_133 <- pc_transcriptomic_df[,1:115] 

pc_transcriptomic_df_133 <- cbind(metadata$status, pc_transcriptomic_df_133)
colnames(pc_transcriptomic_df_133)[1] <- "status"



pc_transcriptomic_df_133$status <- as.factor(pc_transcriptomic_df_133$status) 
s <- sample(nrow(pc_transcriptomic_df_133), nrow(pc_transcriptomic_df_133)*0.7)
training_set <- pc_transcriptomic_df_133[s,]
testing_set <- pc_transcriptomic_df_133[-s,]

model <- glm(status~., training_set, family="binomial", maxit=100)

res <- predict(model, testing_set[,-1], type="response")

table(actual=testing_set$status, predicted=res>10^-1)



