Now we are going to create a heatmap based on the sample counts of the transcripts to identify 1) how well the samples cluster together and 2) the variance of the top 100 genes. 


``` r
library("pheatmap")
library("RColorBrewer")

# rerun dds comparing depth and accounting for locational difference
dds <- DESeqDataSetFromMatrix(countData = gene_abundance_matrix_rounded,
                              colData = metadata,
                              design = ~ location + notes)



# DESeq analysiss
dds <- DESeq(dds)
# save results to res
res <- results(dds)
# summarize results 
summary(res)

# transform scale
vsd <- vst(dds, blind = T)

# plot PCA
plotPCA(vsd, intgroup=c("notes"), pcsToUse = 1:2)
plotPCA(vsd, intgroup=c("notes"), pcsToUse = 2:3)

# get the top most abundant transcripts
select <- order(rowMeans(counts(dds,normalized=FALSE)), decreasing=T)[1:10000]

# save the columns you want to compare
df <- as.data.frame(colData(dds)[,c("location","notes")])

# plot the heatmap
p <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

png("results/2025-04-22-deseq-heatmap-10000-most-abundant.png", width=8.5, height=11, units = "in", res = 300)
print(p)
dev.off()
```

![2025-04-22-deseq-heatmap-10000-most-abundant](https://github.com/user-attachments/assets/bd623ca5-f08c-420d-a57d-64e7b8332f99)


But, what if we only plot the top 100 genes that have the most variance from the mean?

``` r
# transform scale
vsd <- vst(dds, blind = F)

# order by variance and get the top 100 genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)

# extract the information from the vsd file
mat  <- assay(vsd)[ topVarGenes, ]

# difference from the total mean
mat  <- mat - rowMeans(mat)

# get annotation information 
anno <- as.data.frame(colData(vsd)[, c("location","notes")])

# plot the data
p <- pheatmap(mat, annotation_col = anno)

png("results/2025-04-22-deseq-heatmap-100-most-variant.png", width=8.5, height=11, units = "in", res = 300)
print(p)
dev.off()
```

![2025-04-22-deseq-heatmap-100-most-variant](https://github.com/user-attachments/assets/b7cdc0e1-f0fd-48fa-8e12-fe6993654afc)
