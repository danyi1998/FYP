trans_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")

trans_df <- t(trans_df)

trans_df <- as.data.frame(trans_df)

trans_df <- trans_df[vapply(trans_df, function(x) length(unique(x)) > 1, logical(1L))] 



meta <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)

meta <- t(meta)

meta <- as.data.frame(meta)

meta$status <- gsub("bipolar disorder diagnosis: ", "", meta$status)

meta$status <- gsub("BP1", "Patient", meta$status)

meta$status <- gsub("BP2", "Patient", meta$status)

trans_status_df <- cbind(meta$status, trans_df)

colnames(trans_status_df)[1] <- "status"

write.csv(trans_status_df,"C:/Users/65834/Documents/FYP/trans_for_ftest.csv", row.names = FALSE) 
