## While we're getting started... 

Please fill out another round of Minute cards: https://forms.gle/fK2FGG1uUSoaZTSo6

## PICRUST2 Analysis
PICRUST2 is a method to predict functions using phylogenetic placement of ASV sequences - it can help you explore the **functional potential** of a microbial community, which can then (ideally) be tested/validated through other types of -Omics sequencing such as metagenomics or metatranscriptomics. Note: Some people are skeptical of PICRUSt2 results - as they are highly dependent on reference genome databases - but we still think it is a useful (and popular) tool. 
* Github wiki with tutorials - https://github.com/picrust/picrust2/wiki
* Paper to cite: Douglas, G.M., Maffei, V.J., Zaneveld, J.R. et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38, 685–688 (2020). https://doi.org/10.1038/s41587-020-0548-6

![PICRUSt2_flowchart](https://github.com/user-attachments/assets/e072d959-7c83-48e1-848c-607646cc45da)
(Figure from PICRUSt2 wiki)

---

## Predicting Functions Using the PICRUST Workflow

There are several steps to the PICRUST workflow: 

1. Align ASVs to Reference Sequences
2. Place ASVs to Reference Trees
3. Infer Copy Number
4. Metagenome Predications Based on Phylogenetic Placement
5. Infer Pathway Abundances

However, we can easily run PICRUST2 using a single command that wraps all of our steps into one workflow. 

You can only run PICRUST on 16S rRNA metabarcoding datasets, so to show you how to run this pipeline we are going to copy over a projects folder from the instructor data folder. Lets copy the script from the `instructor` data directory onto our personal folder. Remember to replace any instances of `userid` with your account ID. 

```
cp -r /work/mars8180/instructor_data/metabarcoding-datasets/memb-project home/userid/
```

Now we can nano into the scripts and edit our paths. First, lets nano into the script names `07-extract-bio-table-ref-seq.sh`. This script will allow us to export the data (representative sequences and biom table from our QIIME2 artifact files. 

```
cd /home/userid/memb-project/scripts/
nano 07-extract-bio-table-ref-seq.sh
```

```
#!/bin/sh
#SBATCH --job-name="extract-data"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --mail-user=userid@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 07-extract-data.err-%N
#SBATCH -o 07-extract-data.out-%N

module load QIIME2/2024.10-amplicon

REFSEQ=/home/userid/memb-project/results/dada2-rep-seqs.qza
REFTABLE=/home/userid/memb-project/results/dada2-table.qza
SEQOUT=/home/userid/memb-project/results/repseq/
TABLEOUT=/home/userid/memb-project/results/biomtable/

module load QIIME2

qiime tools export \
  --input-path ${REFSEQ} \
  --output-path ${SEQOUT}

qiime tools export \
  --input-path ${REFTABLE} \
  --output-path ${TABLEOUT}
```

Make sure that your paths are pointing to the correct locations. The variable `REFSEQ` and the `REFTABLE` should be path to the ASV representative sequences and the feature table output by the DADA2 script. 

Now we can submit the job to the cluster. This will generate two folders with our exported sequences (`07-repseq/`) and the biom table (`07-biomtable/`).

```
sbatch 07-extract-bio-table-ref-seq.sh
```

I've aleady tried running PICRUST before and I ended up with an error claiming that none of my sequences could be placed to the phylogenetic tree. Out of curiosity, I BLASTed a random sequence and realized that the sequences are in the correct orientation. So, we will use the package `seqtk` to reverse-complement all of our sequences

```
nano 07-rev_complement.sh
```

```
#!/bin/sh
#SBATCH --job-name="rev-com"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --mail-user=userid@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 07-rev-comp.err-%N
#SBATCH -o 07-rev-comp.out-%N

module load seqtk

IN=/home/userid/memb-project/results/07-repseq/dna-sequences.fasta
OUT=/home/userid/memb-project/results/07-biomtable/dna-sequences-rc.fasta

seqtk seq -r ${IN}/ > ${OUT}
```

After you edit your script, you can submit the job. 

```
sbatch 07-rev_complement.sh
```

Now we can use these files to run the picrust2 script and predict functions based on our phylogenetic placement of the ASVs. This time for this scripts, we are going to use 12 threads and 10Gb per thread.

```
nano 08-picrust2.sh
```

```
#!/bin/bash
#SBATCH --job-name="08-picrust2"
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=userid@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e 08-picrust2.err-%N
#SBATCH -o 08-picrust2.out-%N

#Path Variables
SEQ=/home/userid/memb-project/results/07-repseq/dna-sequences-rc.fasta
BIOM=/home/userid/memb-project/results/07-biomtable/feature-table.biom
OUT=/home/userid/memb-project/results/08-picrust/

# Activate Picrust2
module load PICRUSt2

# Picrust2 commands
picrust2_pipeline.py \
  -s ${SEQ} \
  -i ${BIOM} \
  -o ${OUT} -p 12 \
  --stratified \
  --in_traits EC,KO
```

The key output files are:

1. **`EC_metagenome_out`** - Folder containing unstratified EC number metagenome predictions (pred_metagenome_unstrat.tsv.gz), sequence table normalized by predicted 16S copy number abundances (seqtab_norm.tsv.gz), and the per-sample NSTI values weighted by the abundance of each ASV (weighted_nsti.tsv.gz).
2. **`KO_metagenome_out`** - As EC_metagenome_out above, but for KO metagenomes.
3. **`pathways_out`** - Folder containing predicted pathway abundances and coverages per-sample, based on predicted EC number abundances.

## Importing PICRUST files into a phyloseq object 

Phyloseq is agnostic to data types - all it requires is the following

1. Table with the number each function/sequence is observed in each sample
2. Table with a descrption for each row (function/sequence)
3. Table with a description for each sample

Before we can import these files, we need to download these to our results and metadata folder. Typically, we would create a **new** Rproject, but for the sake of time we will save these to our current DDT-project folder. 

Here are the list of files you will have to download from the cluster: 

**Sample Metadata**: `/work/mars8180/instructor_data/metabarcoding-datasets/memb-project/metadata/16S-memb-metadata.txt`
**Function Metadata**: `/work/mars8180/instructor_data/metabarcoding-datasets/memb-project/metadata/ec_traits_pathway.txt`
**Function Biom Table**: `/work/mars8180/instructor_data/metabarcoding-datasets/memb-project/results/08-picrust/KO_metagenome_out/`


After you have downloaded all these files your projects folder, we can import them and convert them into a phyloseq object. 

First, let's use the function `read.delim2` to import our files. 

```
metadata_ec <- read.delim2("metadata/ec_name.txt", sep = "\t", header = F, row.names = 1)
metadata_sample <- read.delim2("metadata/16S-memb-metadata.txt", sep = "\t", header = T, row.names = 1)
ec_predicted <- read.delim2("results/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = T, row.names = 1)
```

Now we need to do some lite data wrangling - we need to make sure that the sample names are consistent across files. We'll use `gsub` to rename some of the colnames from the file `ec_predicted`

```
colnames(ec_predicted) <- gsub("X16Snem", "MEMB.nem.", colnames(ec_predicted))
colnames(ec_predicted) <- gsub("X16SNegCtrl", "Neg.ctrl", colnames(ec_predicted))
colnames(ec_predicted) <- gsub("X16SZymoStandctrl", "Zymo.stand.ctrl", colnames(ec_predicted))
```

Next, we will remove the pattern **EC:** from the file `ec_predicted` and convert the cells to numeric datatype. 

```
rownames(ec_predicted) <- gsub("EC:", "", rownames(ec_predicted))
ec_predicted[,1:304] <- as.numeric(unlist(ec_predicted[,1:304]))
```

This next part will look familiar! We can not convert our files into a phyloseq object. 

```
ec_phy <- otu_table(ec_predicted, taxa_are_rows = TRUE) # notes taxa (AVSs) are as rows
tax_phy <- tax_table(as.matrix(metadata_ec))
samples <- sample_data(metadata_sample)

picrust_phyloseq <- phyloseq(ec_phy, tax_phy, samples)
```

## Analyzing the predicted functional diversity

Now, we can use similar workflows that we have gone through before to analyze our predicted metagenomic data. 

Let's normalize our data by relative abundance. 

```
picrust_phyloseq_relab <- transform_sample_counts(picrust_phyloseq, 
                                                  function(x) x /sum(x) )
```

Now, we will remove samples without any predicted functions (the sum of each column is equal to 0) and the blanks. Typically, we would want to either use decontam or other methods to decontaminate our sample, but for the sake of time we will remove them from the dataset. 

```
picrust_phyloseq_relab_prune <- prune_samples(sample_sums(picrust_phyloseq_relab) >= 1, picrust_phyloseq_relab) # keep samples with at least 1 count
picrust_phyloseq_relab_prune <- subset_samples(picrust_phyloseq_relab_prune, 
                                               FeedingGroup %in% c("1A", "1B", "2A", "2B")) # only keep Feeding Groups 1A, 1B, 2A, 2B 
```

Now, let's use bray-curtis metric to assess dissimilarity of the predicted pathways between samples. 

```
picrust_phyloseq_pcoa <- ordinate(picrust_phyloseq_relab_prune, "PCoA", "bray") 
plot_ordination(picrust_phyloseq_relab_prune, picrust_phyloseq_pcoa, color = "FeedingGroup") + 
  geom_point(size = 3) 
```

If we want to plot the top 20 functions, we can first get a sorted vector and get the top 20 most common functions. Afterwards, let's use this vector to filter our phyloseq object. 

```
top20functions <- names(sort(taxa_sums(picrust_phyloseq_relab_prune), TRUE)[1:20])

#subset phyloseq object to only selected taxa
top20functions_phy <- prune_taxa(family20, picrust_phyloseq_relab_prune) 
```

We will use our favorite R package to plot them as barplots. 

```
plot_bar(top20functions_phy, fill="V2") + 
  facet_grid(~FeedingGroup, space = "free", scales = "free")
```

Finally, we will edit the theme to make the axis text smaller and make the legend more readable. We'll start by making a color palette for our data 

```
library(microViz)
color_list <- distinct_palette(n = 20, pal = "brewerPlus", add = "lightgrey")

plot_bar(top20functions_phy, fill="V2") + 
  geom_bar(aes(color=V2), stat = "identity") +
  facet_grid(~FeedingGroup, space = "free", scales = "free") +
  theme(axis.text.x = element_text(size = 4, vjust = 0.5),
        legend.text = element_text(size=6),
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm')) +
  guides(fill=guide_legend(title="EC Descriptions"),
         color=guide_legend(title="EC Descriptions")) +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list)
```
