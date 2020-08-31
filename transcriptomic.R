transcriptomic_df <- read.table("C:/Users/65834/Documents/FYP/GSE124326_transciptomic.txt", header=TRUE, row.names="gene")
transcriptomic_df <- t(transcriptomic_df)
transcriptomic_df <- as.data.frame(transcriptomic_df)

transcriptomic_df <- transcriptomic_df[vapply(transcriptomic_df, function(x) length(unique(x)) > 1, logical(1L))] 

principal_comps <- prcomp(transcriptomic_df, center = TRUE, scale. = TRUE) 

pc_transcriptomic_df <- principal_comps$x
pc_transcriptomic_df <- as.data.frame(pc_transcriptomic_df)

