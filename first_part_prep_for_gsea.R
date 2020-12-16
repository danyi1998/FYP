reduced_transcriptomic_df <- read.csv(file="C:/Users/65834/Documents/FYP/gene_afterFtest.csv", header=TRUE, sep=",", row.names="samples") 

reduced_transcriptomic_df <- scale(reduced_transcriptomic_df)

metadata <- read.table("C:/Users/65834/Documents/FYP/GSE124326_metadata.txt", header=TRUE, row.names=1)

metadata <- t(metadata)

metadata <- metadata[!(row.names(metadata) %in% c("X12R1998", "X13R1190")), ]

reduced_transcriptomic_df <- reduced_transcriptomic_df[!(row.names(reduced_transcriptomic_df) %in% c("X12R1998.counts", "X13R1190.counts")), ]

metadata <- as.data.frame(metadata)
reduced_transcriptomic_df <- as.data.frame(reduced_transcriptomic_df)

reduced_transcriptomic_df <- cbind(metadata[,1:2], reduced_transcriptomic_df)

reduced_transcriptomic_df <- reduced_transcriptomic_df[,-2]

reduced_transcriptomic_df <- reduced_transcriptomic_df[reduced_transcriptomic_df$status != "bipolar disorder diagnosis: Control", ]

reduced_transcriptomic_df <- reduced_transcriptomic_df[,-1]

reduced_transcriptomic_df <- t(reduced_transcriptomic_df)

write.csv(reduced_transcriptomic_df,"C:/Users/65834/Documents/FYP/first_part_prep_for_gsea.csv", row.names = TRUE) 

