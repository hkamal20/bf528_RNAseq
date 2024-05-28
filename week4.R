#install.packages("DESeq2")
#install.packages("tidyverse")

library(DESeq2)
library(tidyverse)
library(dplyr)

counts <- as.matrix(read.csv('results/all_non_zero_verse_counts.csv', row.names = 'gene'))
#head(counts)

coldata <- data.frame(samples = colnames(counts), time = c(rep('AD', 2), rep('P0', 2), rep('P4', 2), rep('P7', 2)), row.names = 'samples')
coldata$time <- as.factor(coldata$time)
coldata$time <- relevel(coldata$time, ref='P0')
#coldata


dds <- DESeqDataSetFromMatrix(countData = counts, colData=coldata, design = ~time)
dds <- DESeq(dds)
res <- results(dds, contrast = c('time', 'AD', 'P0'))
results_df <- as.data.frame(res)
#head(results_df)
gene_names <- rownames(results_df)
results_df <- data.frame(gene_id = gene_names, results_df)
#head(results_df)
results_tib <- as_tibble(results_df)
#head(results_tib)


id2 <- read.csv('results/id2gene.txt', sep = '\t')
colnames(id2) <- c('gene_id', 'gene_name')

#head(id2) 


merged_df <- merge(results_tib, id2, by.x = 'gene_id', by.y = 'gene_id', all.x = TRUE)
filt_df <- filter(merged_df, !is.na(log2FoldChange)) 
sort <- filt_df[order(-filt_df$log2FoldChange), ]
rank_df <- setNames(sort$log2FoldChange, sort$gene_name)

head(rank_df)
head(filt_df)
#resOrdered <- res[order(res$pvalue),]
#resOrdered %>% as_tibble(rownames = 'geneid')

#res %>% as_tibble(rownames = 'geneid') %>% left_join(id2, by='geneid') %>% arrange(padj) %>% select(geneids, genenames, padj, log2FoldChange)