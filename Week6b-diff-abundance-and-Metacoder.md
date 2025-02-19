### Week 6b Other Ecological Statisics and R Packages

## About our samples - study design and key questions:

Previous studes on DDT barrel dumpsites:
* Neira et al. 2024 Waste barrel contamination and macrobenthic communities in the San Pedro Basin DDT dumpsite, *Marine Pollution Bulletin*, 203: 116463 - https://www.sciencedirect.com/science/article/pii/S0025326X24004405 (macrofauna patterns & impacts)
* Kivenson et al. 2019 Ocean dumping of containerized DDT waste was a sloppy process, *Environmental Science and Technology*, 53(6):2791-2980 - https://pubs.acs.org/doi/10.1021/acs.est.8b05859 (detailed background information about DDT barrel dumpsites + what we know science-wise)

---

<img width="1013" alt="Screenshot 2025-02-18 at 8 39 32 AM" src="https://github.com/user-attachments/assets/063a6ca1-e142-4d85-af59-bafb3695b0b6" />

---

<img width="2531" alt="Screenshot 2025-02-18 at 8 40 06 AM" src="https://github.com/user-attachments/assets/9cdb1343-2564-409f-a730-095ac4691ecd" />

---

The Bik Lab DDT work is focusing on:
* Assessing eukaryotic biodiversity patterns using an 18S rRNA metabarcoding approach (using DNA extraction aliquots we obtained from the Jensen lab at UCSD who are doing parallell 16S rRNA analysis of bacterial/archael communities)
* Looking for nematode "bioindicator" taxa indicative of DDT pollution (e.g. taxa prevalent on/inside the white chemical ring)

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

First, lets focus on just the first 2cm sediment so that we can better compare the different Sites.

```
barrel_one_0_2 <- subset_samples(phylo_obj_tree_sans_contam, Core_Fraction=="0_2") # subset our data 
barrel_one_0_2 <- prune_taxa(taxa_sums(barrel_one_0_2) > 0, barrel_one_0_2) # remove ASVs not present
```

Now, lets subset our dataset so we only have ASVs Identified as Nematoda and then agglomerate that taxa at genus level 

```
nematoda_phy <- subset_taxa(barrel_one_0_2, Taxon14=="D_13__Nematoda")
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
aldex_taxa_info <- aldex_taxa_info %>% rownames_to_column(var = "OTU") # change rownames to column named OTU
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

Now, let's get the OTU table from our phyloseq object and merge it with out ALDEX2 results

```
gen_otu_table_location <- data.frame(phyloseq::otu_table(nematoda_phy_taxon22)) # extract otu table
gen_otu_table_location <- rownames_to_column(gen_otu_table_location, var = "OTU") # covert rownames to column

sig_aldex2_gen_result_location_tax <- left_join(sig_aldex2_gen_result_location, gen_otu_table_location) # combine the tables
```

Let's change the rownames from the OTU hash to the assigned Genus

```
rownames(sig_aldex2_gen_result_location_tax) <- sig_aldex2_gen_result_location_tax$Taxon22 # change Taxon22 (Genus) to rownames
sig_aldex2_gen_result_location_tax <- sig_aldex2_gen_result_location_tax[, -(1:27)] # remove tax info from otu table
```

Now, lets log transform the otu table and convert them the reads to a z-score.

```
shsk_gen_czm_location <- (apply(shsk_gen_czm_location, 1, function(x){log(x+1) - mean(log(x+1))})) # log transform
Z.Score.gen_location <- scale(t(shsk_gen_czm_location)) # scale the dataset
```

Now let's prepare for our visualization by order our dataset and creating our heatmap annotations and labels

```
df <- sample_data(nematoda_phy_taxon22) # extract metadata from phyloseq object
list=df[order(df$Site, decreasing=T),] # order based on Site so we can group them together
sample=rownames(list) # extract the rownames from the list (SampleID)
heatmap_annotation_top = HeatmapAnnotation(
  Location = df$Site) # We will annotate by site 
```

Now, we are going to create a scale so we easily visualize our z-score (high z-score is green, low z-score is brown)

```
col_matrix <- brewer.pal(6, "BrBG")
```

Finally, after spending all this time data wrangling, we can use the function `Heatmap` to plot out data. The file `Z.Score.gen_location` has our Z-scores. The paramter `column_order` will allow use to group samples from the same Sites together and `col` lets use set the color palette. 

```
hm_gen_location <- Heatmap(Z.Score.gen_location, name = "Z-score", col = col_matrix,
                           column_order = sample, 
                           top_annotation = heatmap_annotation_top,
                           cluster_columns = F, 
                           column_names_gp = gpar(fontsize = 6))
```

How many nematode genera are differentially abundant in our dataset?

## Visualizing taxonomy with Metacoder

Metacoder is a package that allows us to visualize taxonomic (or gene ontology) data using "heat trees, instead of the typical barplots. It can be a great way to get a different kind of visual overview of your dataset, and identify taxa/genes that warrant further exploration 
* GitHub repo (documentation and tutorials): https://grunwaldlab.github.io/metacoder_documentation/
* Publication to cite: Foster ZSL, Sharpton TJ, Grünwald NJ (2017) Metacoder: An R package for visualization and manipulation of community taxonomic diversity data. PLoS Comput Biol 13(2): e1005404. https://doi.org/10.1371/journal.pcbi.1005404

---

<img width="1032" alt="Screenshot 2025-02-19 at 6 30 20 PM" src="https://github.com/user-attachments/assets/1cb82493-9c94-4016-b05f-faaabdf4fd12" />

---

First, we need to install this package using the following command. 

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
