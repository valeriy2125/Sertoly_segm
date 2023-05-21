install.packages(DESeq2)
library("DESeq2")

# Read the data with raw counts from the table.
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
condition <- c('1_3', '1_3', '4_6', '1_3', '4_6', '4_6') #vector of column names for the data frame

# design table
colData <- data.frame(row.names=colnames(countData), condition=factor(condition, levels=c('1_3','4_6')))
colData

dataset <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~condition)

# Deseq dataset creation
dds <- DESeq(dataset)

result <- results(dds, contrast=c('condition','4_6','1_3'))
result <- result[complete.cases(result),]  #remove any rows with NA
head(result)

write.csv(result, "result_old.csv")


