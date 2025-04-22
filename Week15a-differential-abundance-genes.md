## Differential Gene Abundance Analysis
We are going to use the DESeq2 R pacakge [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8]) to conduct a differential gene abundance analysis to identify the genes that are more abundant in the surface water compared to the deep chlorophyll maximum. 

First, we need to download three files from the computing cluster.

1. The salmon merged quant files
2. The emapper gene annotations
3. The metadata file

We can download them to our computer using the following commands:

```
scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metatranscriptome-datasets/results/13-eggnog/tara.emapper.annotations-contig-gene.txt .

scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metatranscriptome-datasets/results/10-salmon/salmon-all-sample-quant-counts.txt .

scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metatranscriptome-datasets/metadata/2025-04-05-metatranscriptomics-tara-oceans-metadata.txt
```

```r
#### load your packages 
#### install them if not previosly installed
library(tximport)
library(DESeq2)
library(tidyverse)
library(data.table)


#### METADATA FILE #####
# import the metadata file
metadata <- read.delim2("metadata/2025-04-05-metatranscriptomics-tara-oceans-metadata.txt", row.names = 1)
#rownames(metadata) <- metadata$sampleID
metadata <- metadata[ order(row.names(metadata)), ] # order the metadata file based on the row name



#### SALMON FILE ##### 
# import salmon quants
salmon_data <- fread("data/salmon-all-sample-quant-counts.txt", header = T)
# coerce to a dataframe
setDF(salmon_data) 
# change the contig name to row name
rownames(salmon_data) <- salmon_data$Name
# remove the column 
salmon_data$Name <- NULL 



#### EMAPPER FILE #####
# import the emapper annotations
tx2gene <- fread("data/tara.emapper.annotations-contig-gene.txt", header = F, skip = 5)
# coerce to a dataframe
setDF(tx2gene)
# remove ".p#" after contig
tx2gene$V1 <- gsub(".p[0-9]", "", tx2gene$V1)




#### MERGE DATA ####
# merge the data by contig
merged_data <- merge(salmon_data, tx2gene, by.x = "row.names", by.y = "V1")
# remove row.names column
merged_data$Row.names <- NULL




#### EXTRACT SAMPLE NAMES ####
sample_cols <- names(salmon_data)




#### AGGREGATE DATA ####
# aggregate data by gene 
gene_abundance_dt <- aggregate(. ~ V2, data = merged_data[, c(sample_cols, "V2")], FUN = sum)
# make gene name the row name
rownames(gene_abundance_dt) <- gene_abundance_dt$V2
# remove the column V2
gene_abundance_dt$V2 <- NULL
# convert to a matrix 
gene_abundance_matrix <- as.matrix(gene_abundance_dt)
# round the matrix 
gene_abundance_matrix_rounded <- round(gene_abundance_matrix)
# convert to DESeq Data Set 
dds <- DESeqDataSetFromMatrix(countData = gene_abundance_matrix_rounded,
                              colData = metadata,
                              design = ~ notes)




#### DESeq analysis #####
dds <- DESeq(dds)
resultsNames(dds)
# save results to res
res <- results(dds)
# summarize results 
summary(res)

#### VISUALIZE ####
# let's visualize using a Volcano plot
p <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     title = 'Volcano plot',
                     pCutoff = 0.05,
                     FCcutoff = 2,
                     pointSize = 1,
                     labSize = 2)
p

#### SAVE ####
# save as png
png("results/2025-04-22-salmon-emapper-volcano-plots.png", width=11, height=8.5, units = "in", res = 300)
print(p)
dev.off()
```
![2025-04-22-salmon-emapper-volcano-plots](https://github.com/user-attachments/assets/46e08793-2040-44d1-b90b-f01a96213525)

