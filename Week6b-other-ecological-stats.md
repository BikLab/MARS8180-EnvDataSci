### Week 6b Other Ecological Statisics and R Packages

Cover PICRUST + other tools this week (we had student requests for this)

## Calculating Differentially Abundant Nematode Genera

Aldex2 was first developed to analyzed RNA-seq datasets; however, it has been shown to work well with compositional datasets, including metabarcoding data. We will use this package to identify nematode genera that are differentially abundance in sample areas. 

First, we need to install a few packages

```
BiocManager::install("ALDEx2")
BiocManager::install("RColorBrewer")
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")
```

After, lets load them into our R environment 

```
library(ALDEx2)
library(tidyr)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
```

Now, lets subset our dataset so we only have ASVs Identified as Nematoda and then agglomerate that taxa at genus level 

```
nematoda_phy <- subset_taxa(phylo_obj_tree_sans_contam, Taxon14=="D_13__Nematoda")
nematoda_phy_taxon22 <- tax_glom(nematoda_phy, taxrank="Taxon22")
```

Let's run the aldex2 pipeline using a kruskall-wallis test 

```
aldex2_sample_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(nematoda_phy_taxon22)),
                                     phyloseq::sample_data(nematoda_phy_taxon22)$Site_Area,
                                     test="kw")
```

We will have to convert our `otu_table` into a dataframe and make the rownames into a column named **OTU**.

```
gen_otu_table_location <- data.frame(phyloseq::otu_table(nematoda_phy_taxon22)) # create dataframe
gen_otu_table_location <- rownames_to_column(gen_otu_table_location, var = "OTU") # convert rownames to column
```

Let's do the same thing for our Aldex output

```
aldex2_sample_site$OTU <- row.names(aldex2_sample_site)
row.names(aldex2_sample_site) <- NULL
```

We are going to convert the taxonomy to a dataframe

```
aldex_taxa_info <- data.frame(tax_table(nematoda_phy_taxon22)) # convert taxonomy table to dataframe
aldex_taxa_info <- aldex_taxa_info %>% rownames_to_column(var = "OTU") change rownames to column named OTU
```

We will filter our results, to only keep the ASV/OTUs with a significant kw.ep and append taxonomy information to it.

```
# Filter aldex2 results by sig kw.ep and join the taxanomic information
aldex2_sample_site_sig <- aldex2_sample_site %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
  
sig_aldex2_gen_result_location <- left_join(aldex2_sample_site_sig, aldex_taxa_info) # append taxonomy info
```
