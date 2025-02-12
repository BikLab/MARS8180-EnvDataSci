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

Let's convert our Aldex output into a dataframe and make the rownames into a column 

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

How many nematode genera are differentially abundant in our dataset?

## Visualizing taxonomy with Metacoder

Metacoder is a package that allows us to visualize taxonomic data using heat trees, instead of the typical barplots. First, we need to install this package using the following command. 

```
BiocManager::install("metacoder")
library(metacoder)
```

Now, we need to transform our counts to relative abundance. 

```
phylo_obj_tree_sans_contam_ra <- transform_sample_counts(phylo_obj_tree_sans_contam, function(x) x / sum(x))
```

Now, we want to focus on comparing two samples 1) the 2cm core fraction at the North Barrel Site and 2) the 2cm core fraction at the San Diego Trough. We can do this by subsetting our dataset and then removing ASVs that are not present in this subsetted dataset. Let's start with the 2cm core fraction at the North Barrel Site. 


```
barrel_one_0_2 <- subset_samples(phylo_obj_tree_sans_contam, Site=="DDT_Barrel_Site_1_North" & Core_Fraction=="0_2") # subset our data 
barrel_one_0_2 <- prune_taxa(taxa_sums(barrel_one_0_2) > 0, barrel_one_0_2) # remove ASVs not present
```

Let's subset the dataframe so we are only comparing the ASVs assigned to **Nematoda**

```
barrel_one_nematoda <- subset_taxa(barrel_one_0_2, Taxon14=="D_13__Nematoda") # subset nematoda
tax_table(barrel_one_nematoda) <- tax_table(barrel_one_nematoda)[,c(14, 15:17,21,22)] # only keep major taxonomic heirarchies starting at the phylum level. 
```

Let's convert the phyloseq object into a specific format called a `taxmap`

```
b1_obj <- parse_phyloseq(barrel_one_nematoda)
```

Finally, we can plot our results

```
b1_obj_metacoder <- heat_tree(b1_obj,
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs)
```

#### Lets do the same thing for the second sample set (2cm core fraction at the San Diego Trough)

```
sd_trough <- subset_samples(phylo_obj_tree_sans_contam, Site=="San_Diego_Trough" & Core_Fraction=="0_2") # subset samples 
sd_trough <- prune_taxa(taxa_sums(sd_trough) > 0, sd_trough) # remove ASVs not present

sd_trough_nematoda <- subset_taxa(sd_trough, Taxon14=="D_13__Nematoda") # only keep Nematoda ASVs
tax_table(sd_trough_nematoda) <- tax_table(sd_trough_nematoda)[,c(14, 15:17,21,22)] # only keep major taxonomic heirarchies

sd_obj <- parse_phyloseq(sd_trough_nematoda) # convert the phyloseq object to taxmap

sd_obj_metacoder <- heat_tree(sd_obj,
                              node_label = taxon_names,
                              node_size = n_obs,
                              node_color = n_obs) # visualize the results

```
