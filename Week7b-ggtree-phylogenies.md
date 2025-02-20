### Week 7b Phylogenetic Visualizations using ggtree

One type of underutilized visualization in -Omics dataset is phylogenetic trees - you can easily build these with rRNA data (18S for eukaryotes / 16S for bacteria and archaea), and phylogenies will give you a lot of additional information about the evolutionary history of microbial assemblages.

Our favorite package is ggtree, and the supplemental package ggtree extra, software papers as follows:

* G Yu, DK Smith, H Zhu, Y Guan, TTY Lam*. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution. 2017, 8(1):28-36. - https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628

* Xu S, Dai Z, Guo P, Fu X, Liu S, Zhou L, Tang W, Feng T, Chen M, Zhan L, Wu T, Hu E, Jiang Y, Bo X, Yu G (2021). “ggtreeExtra: Compact visualization of richly annotated phylogenetic data.” Molecular Biology and Evolution, 38, 4039-4042. ISSN 0737-4038, doi:10.1093/molbev/msab166, https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msab166/6294410.

An example ggtree figure from a recent Bik Lab paper focused on nematode microbiomes (phylogeny of host nematode species, and bar plots of key bacterial microbiome taxa):

![Figure_5](https://github.com/user-attachments/assets/4b3ca216-75ef-4e83-b00d-c67af6dd0afb)

(Figure 5 from Pereira TJ, De Santiago A, Bik HM (2024). Soil properties predict below‐ground community structure, but not nematode microbiome patterns in semi‐arid habitats. Molecular Ecology, 33(18), e17501.- https://onlinelibrary.wiley.com/doi/full/10.1111/mec.17501

---

## Utilizing ggtree to analyze metabarcoding data
There are several software packages that extend the use of ggplot2 to visualize -omics datasets. Lets use ggtree to visualize the abundance of different nematode genera across Sample Sites. 

As always, lets first install the necessary packages:

```
install.packages("ggstar")
install.packages("ggnewscale")
install.packages("dplyr")
install.packages("ggtree")
install.packages("ggtreeExtra")
install.packages("ggplot2")

library(ggstar)
library(ggnewscale)
library(dplyr)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)

```

Next, we are going to subset are data and agglomerate at the genus level to focus solely on marine nematodes. We are going to start with our phyloseq object `phylo_obj_tree_sans_contam`.

```
nematoda_phy <- subset_taxa(phylo_obj_tree_sans_contam, Taxon14=="D_13__Nematoda") # subset nematoda
nematoda_phy_taxon22 <- tax_glom(nematoda_phy, taxrank = "Taxon22") # agglomerate genus (Taxon22)
```

Afterwards, we can merge the samples based on the **Site** and normalize our reads to relative abundance.

```
nematoda_phy_merge_samples <- merge_samples(nematoda_phy_taxon22, "Site", fun = mean)
nematoda_phy_merge_samples_ra <- transform_sample_counts(nematoda_phy_merge_samples, function(x) x / sum(x) * 100)
``` 

We are going to melt both of those phyloseq object:

```
melt_sample_reads <- as.data.frame(psmelt(nematoda_phy_merge_samples)) # "melt" data and force to a dataframe
melt_sample_ab <- as.data.frame(psmelt(nematoda_phy_merge_samples_ra)) # "melt" data and force to a dataframe
```

Let's subset our data so we just have the information we want to plot. This part is important - sometime ggtreeExtra does not play well if you are integrating multiple dataframes with the same columns. To avoid this, we aare going to rename abundance to relative abundance for the normalized data. For the non-normalized data, we will rename Abundance to Reads. 

```
melt_sample_ab <- melt_sample_ab %>% select(OTU, Relab = Abundance, Site=Sample) # keep three columns and rename two
melt_sample_reads <- melt_sample_reads %>% select(OTU, Reads = Abundance, Sample, Taxon19, Taxon16) # keep 5 columns and rename 1
```

Now, we need to extract our phylogenetic tree from our data

```
tree <- phy_tree(nematoda_phy_taxon22)
```

Let's build our rough circularized phylogenetic tree

```
p <- ggtree(tree, size=0.1, open.angle=5, branch.length = "none", layout = "fan") 
p
```

**What are the paramters/flags doing?**

We are going to use the `melt_sample_reads` dataframe to color and shape the tree tips according to taxa. Additionally, there size will correlate to the total number of reads in the dataset. We will use a pipe function `%<+%` that will allow use to use the file as the input. 

```
p2 <- p %<+% melt_sample_reads + # use this file as the input 
  geom_tippoint(mapping = aes(color = Taxon19, shape = Taxon16, size = Reads))  + # set variable for color, shape, and tip size
  scale_color_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080", # change colors
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970", "black")) 
p2 
```

Now, lets add more information - lets include a heatmap that correspond to the relative abundance in each Sample Site. To do this, we will use the file `melt_sample_ab`. 

```

p3 <- p2 + new_scale_fill() + new_scale_color() + # reset the color and fill scales
  geom_fruit(data = melt_sample_ab, # input file
             geom=geom_tile, # type of graph 
             mapping=aes(y=OTU, x=Site, alpha=Relab, fill=Site), # input the coordinates (x and y variables) and the color and transparancy
             offset = 0.02, size = 0.05, pwidth = 0.1, color="grey") +
  scale_alpha_continuous(range=c(0, 1), limits=c(0, 100), breaks=seq(0,100, by=10)) + # set the ranges for the fill
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000", # change colors
                             "#800000", "#006400","#800080","#696969")) +
  theme(legend.position=c(1.15, 0.5), 
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.spacing.y = unit(0.02, "mm"),
        legend.key.height= unit(2, 'mm'),
        legend.key.width= unit(4, 'mm')) 
  
p3
```





