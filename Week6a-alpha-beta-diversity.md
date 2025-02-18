### Week 6a Alpha and Beta Diversity Analyses

#### Statistical analyses of metabarcoding data

Some key points from Amy Willis's lectures at the STAMPS 2024 course at MBL: https://github.com/mblstamps/stamps2024/wiki#17

---

<img width="814" alt="Screenshot 2025-02-02 at 2 07 34 PM" src="https://github.com/user-attachments/assets/b737bb83-8771-4ee4-80e0-f1aade176d43" />
<img width="813" alt="Screenshot 2025-02-02 at 2 07 21 PM" src="https://github.com/user-attachments/assets/a8e990a8-3c84-4a68-9b46-973457236b5e" />

---

<img width="822" alt="Screenshot 2025-02-02 at 2 07 58 PM" src="https://github.com/user-attachments/assets/e7ee22c3-5c0c-405f-9de6-25b3421c7231" />

---

Mbareche et al. 2020 - https://www.mdpi.com/2075-1729/10/9/185

<img width="963" alt="Screenshot 2025-02-11 at 10 09 57 AM" src="https://github.com/user-attachments/assets/b2a0994f-cab7-44a6-8797-f9ec98d62991" />

---

Gloor et al. (2017) Microbiome datasets are compositional: And this is not optional, *Frontiers in Microbiology*, 8:2224 - https://doi.org/10.3389/fmicb.2017.02224 - **you must read this paper (especially before the midterm)**

<img width="542" alt="Screenshot 2025-02-11 at 10 18 06 AM" src="https://github.com/user-attachments/assets/836743f1-627f-4b27-86ea-9ecb5671a849" />

---

## Generating taxonomy barplots

Taxonomy barplots are a common way to visualize our community composition. We can do this in phyloseq, which wraps ggplot2, to create these visualizations. That means we can use the ggplot2 language to further enchange our graphs. 

We will focus on visualizing the community structure of the marine nematodes since we have well-curated strings. So first, we need to subset our dataframe to only include nematodes and then transform/normalize our counts to relative abundance. 

```
phylo_obj_tree_sans_contam_nematoda <- subset_taxa(phylo_obj_tree_sans_contam, Taxon14=="D_13__Nematoda") # subset for marine nematodes
phylo_obj_tree_sans_contam_nematoda_ra <- transform_sample_counts(phylo_obj_tree_sans_contam_nematoda, function(x) x / sum(x)) # convert to relative abundance
```

Now, we can plot using the `plot_bar` command. I will include other ggplot functions to make sure our figure looks good. This includes `facet_grid` which will separate our samples by Site and the flag `free-scales` will not plot samples in the facet if they do not belong to the Site.  If you want to see how the command reacts without these parameters, remove them :). 

```
plot_bar(phylo_obj_tree_sans_contam_nematoda_ra, "Sample", fill="Taxon17") + # Sample is x-axis
  geom_bar(aes(color=Taxon17), stat="identity") + # specifies the type of graph and color
  facet_grid(~Site, scales = "free", space = "free_x") # lets separate our samples by Site
```

Our titles are too long and unreadable, so lets edit them to make this publishable. First, lets create a dictionary and create "long" form labels.

```
long_form_labels <- c(
  DDT_Barrel_Site_1_North = "DDT Barrels North",
  DDT_Barrel_Site_2_South = "DDT Barrels South",
  Patton_Ridge_South = "Patton Ridge",
  `40Mile_Bank` = "40Mile Bank",
  Lasuen_Knoll = "Lasuen Knoll",
  San_Diego_Trough = "San Diego Trough",
  no.data = "controls")
```

Now, we need to match them to each samples and save them to a new column

```
idx <- match(sample_data(phylo_obj_tree_sans_contam_nematoda_ra)$Site, names(long_form_labels)) # match to each Sample
sample_data(phylo_obj_tree_sans_contam_nematoda_ra)$Site_long <- long_form_labels[idx] # create a new column in our phyloseq metadata
```

Finally, we can plot our data so its readable. The `labeller` flag in the facet command will create line breaks. The `theme()` command allows us to change the font sizes and position. 

```
plot_bar(phylo_obj_tree_sans_contam_nematoda_ra, "Sample", fill="Taxon17") +
  geom_bar(aes(color=Taxon17), stat="identity") +
  facet_grid(~Site_long, scales = "free", space = "free_x", labeller = labeller(Site_long = label_wrap_gen(1))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        strip.text.x = element_text(size = 8, angle = 90))
```

## Alpha Diversity, Rarefaction, and Normalization

Key Papers:
* Willis AD (2019) Rarefaction, Alpha Diversity, and Statistics, *Frontiers in Microbiology* -https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.02407/full
* McMurdie, P. J., and Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. PLoS Comput. Biol. 10:e1003531 - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531

Alpha diversity metrics summarize the structure of an ecological community with respect to its richness (number of taxonomic groups), evenness (distribution of abundances of the groups), or both. Because many perturbations to a community affect the alpha diversity of a community, summarizing and comparing community structure via alpha diversity is a ubiquitous approach to analyzing community surveys. In microbial ecology, analyzing the alpha diversity of amplicon sequencing data is a common first approach to assessing differences between environments. (Willis 2019 - see paper link above)

**Rarefaction** - Subsampling your overall dataset, by selecting X number of observations (e.g. ASVs) from each sample, regardless of the overall sequencing effort. For example, if you have two sample sites, Beach A (1000 ASVs and 1 million sequence reads), and Beach B (250 ASVs and 200,000 sequencing reads), you could rarify this dataset by randomly subsampling 100 ASVs from each sample site. Rarefaction is a common paremeter you will see in the scripts/code we run in this class. **Note: Rarefaction can be contentious in eDNA studies, but many people still use it - see above literature refs**

**Normalization** - Using proportional abundances of ASVs, instead of the count of sequence reads for that ASV. To normalize a metabarcoding dataset, you typically divide the raw sequence counts for each ASV in a sample by the total number of sequences in that sample. For example, at our sample site Beach A, ASV_347 has an absolute abundance of 50,000 sequence reads. If Beach A has 1 million sequence reads overall, then the normailized abundance of ASV_347 would be 50,000/1,000,000 --> 0.05 or 5% 

Something to think about (annotated image from Willis et al. 2019): 

<img width="1085" alt="Screenshot 2025-02-13 at 10 15 42 AM" src="https://github.com/user-attachments/assets/4230fd18-b62d-4305-9732-0b0bd758b90f" />

---

## Calculating Alpha Diversity

Key Alpha Diversity metrics and visualizations:
* Taxonomy bar charts
* Species Richness
* Evenness
* Chao1, Shannon, Simpson Diversity Indexes

Now that we have created our phyloseq object, we are going to calculate Observed ASVS - a simple measure of alpha diversity. We can do this by using the `phyloseq::estimate_richness` function. We usually do not need to to specify the package, but this might be useful if you are using packages that have functions with the same name. 

```
alpha_div_observed <- phyloseq::estimate_richness(phylo_obj_tree_sans_contam, measures = "Observed") # calculate alpha diversity
head(alpha_div_observed) # view the first few lines of our R object
```

```
         Observed
DDT.1.1       356
DDT.1.2       384
DDT.10.1      352
DDT.10.2      406
DDT.11.1      473
DDT.11.2      480
```

We were able to calculate diverstiy, but we lost our metadata associated with our samples. However, we can access this metadata using the `sample_data` function.

```
metadata_df <- phyloseq::sample_data(phylo_obj_tree_sans_contam) # access metadata from the phyloseq object
head(metadata_df) 
```

We can put this together to calculate the observed diversity (number of ASVs) and create a dataframe that includes our metadata.

```
alpha_div_observed_metadata <- data.frame( # save as a dataframe
  phyloseq::sample_data(phylo_obj_tree_sans_contam), # access metadata 
  "Observed" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam, measures = "Observed")) # calculate Observed diversity and save it to a column names "Observed"
```  

Now we can use ggplot to plot our results

```
alpha_div_observed_plot <- ggplot(alpha_div_observed_metadata, aes(x=Barrel_ID, y=Observed)) + # what are we plotting
  geom_boxplot() + # make it into a boxplot plot
  theme_minimal() # and lets make sure its "pretty"  
```

We can modify the axis labels to make it more readable. We are going to tilt the labels to 45 degrees, change the font size to 8, and adjust the horizontal justification

```
alpha_div_observed_plot <- ggplot(alpha_div_observed_metadata, aes(x=Site_Area, y=Observed)) + # what are we plotting
  geom_boxplot() + # make it into a boxplot plot
  theme_minimal() + # and lets make sure its "pretty"
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # change x-axis labels
```

## Calculating Several Alpha Diversity Metrics
Lets do the same by calculating alpha diversity using commonly used metrics, such as number of reads, shannon diversity, simposon, and eveness.

```
alpha_div_observed_metadata <- data.frame(
  phyloseq::sample_data(phylo_obj_tree_sans_contam), # get metadata
  "Reads" = phyloseq::sample_sums(phylo_obj_tree_sans_contam), # number of reads
  "Observed" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam, measures = "Observed"), # count observed ASVs
  "Shannon" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam, measures = "Shannon"), # Calculate Shannon Diversity
  "InvSimpson" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam, measures = "InvSimpson")) # calculate InvSimpson
```

Now lets calculate Eveness by dividing Shannon by the Log of the Observed ASVs

```
alpha_div_observed_metadata$Evenness <- alpha_div_observed_metadata$Shannon/log(alpha_div_observed_metadata$Observed)
head(alpha_div_observed_metadata)
```

## Beta Diversity - Background and Resources

Key Beta Diversity metrics and visualizations:
* Bray-Curtis - considers ASV presence/absence AND abundance, gives more consideration to abundant ASVs
* Jaccard - only consideres ASV presence/absence
* Canberra - simliar to Bray-Curtis, but gives more consideration to rare taxa (low abundance ASVs) by treating all ASVs equally
* Unifrac (weighted/unweighted) - take phylogenetic distance into account, can be based on presence/absence of ASVs only (unweighted - recommended), or also additionally incorporate normalized abundances (weighted - less commonly used)
* nMDS plots

If you really want to get into the weeds with some of these metrics, here are some recommended resources:
* Roberts DW (2017) Distance, dissimiliarity, and mean-variance ratios in ordination, Methods in Ecology and Evolution, 8:1398-1407 - https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210X.12739
* QIIME2 online documentation also has some pretty good explanations of different metrics (alpha/beta diversity) - https://docs.qiime2.org/
* See also some of the phyloseq tutorials (with links to key literature) - https://joey711.github.io/phyloseq/

## Ordination Crash Course (class request)

A really good overview of ordination types - https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/types-of-ordination-methods/ (**NOTE:** This whole book is freely available online, and highly recommended as a deep dive - Baker (2024) Applied Multivariate Statistics in R - https://uw.pressbooks.pub/appliedmultivariatestatistics/)

Ordination is a _**"graphical representation of the similiarity of sampling units and/or attribues in resemblance space"**_ (Wild 2010, pg. 35), quote from Baker textbook

All ordination methods are aimed at **identifying hidden patterns** in a dataset and **reducing the complexity of high-dimensional data**. This is often done by carrying out complex mathematical calculations and comparison of samples, and then "squishing" this data down into a 2D plot, by maximizing the usefulness of the visualization (e.g. choosing the aspect of the data which differentiates sample groupings the most, such as a PCoA axis).  

There are many types of ordination methods, including (but not limited to):
* **Principal Component Analysis (PCA)** - "a statistical technique used to reduce the dimensionality of a dataset while retaining most of its variability. It is a linear transformation method that converts the original set of variables into a new set of linearly uncorrelated variables, called principal components (PCs), which are sorted in decreasing order of variance." It primarily uses Euclidian distance - https://r.qcbs.ca/workshop09/book-en/principal-component-analysis.html (**NOTE:** in QIIME2 you choose the distance matrix used to build your PCoAs, such as Bray-Curtis similarity or the Unifrac metric which measures phylogenetic distance betwen samples/species). 
* **Correspondence Analysis (CA)** - "preserves Chi2 distances" and may be a better option for data with "long ecological gradients", good for when underlying assumptions of PCA are violated by the dataset characteristics - https://r.qcbs.ca/workshop09/book-en/correspondence-analysis.html
* **Principal Coordinates Analysis (PCoA)** - an "unconstrained ordination" where "points are added to plane space one at a time using Euclidean distance (or whatever distance (dissimilarity) metric you choose)" - https://r.qcbs.ca/workshop09/book-en/principal-coordinates-analysis.html
* **Nonmetric Multidimentional Scaling (NMDS)** - in NMDS, "the priority is not to preserve the exact distances among sites, but rather to represent as accurately as possible the relationships among objects in a small and number of axes (generally two or three) specified by the user", and so "the biplot produced from NMDS is the better 2D graphical representation of between-objects similarity: dissimilar objects are far apart in the ordination space and similar objects close to one another" - https://r.qcbs.ca/workshop09/book-en/nonmetric-multidimensional-scaling.html
* Other ordination types mentioned in class: **Detrended Correspondence Analysis (DCA)**, **ReDundancy Analysis (RDA)**, **Distance-based Redundancy Analysis (db-RDA)** also known as **Canonical Analysis of Principal Coordinates (CAP)** - see Baker texbook (above link) for great explanations of all of these and more!

*Choosing between PCA and PCoA can be tricky, but generally PCA is used to summarize multivariate data into as few dimensions as possible, whereas PCoA can be used to visualize distances between points. PCoA can be particularly suited for datasets that have more columns than rows. For example, if hundreds of species have been observed over a set of quadrats, then a approach based on a PCoA using Bray-Curtis similarity may be best suited.* (Quote from: https://r.qcbs.ca/workshop09/book-en/principal-coordinates-analysis.html) 

Goodrich et al. (2016) Conducting a Microbiome Study, Cell, 158(2):250-262 - [https://pmc.ncbi.nlm.nih.gov/articles/PMC5074386/](https://doi.org/10.1016/j.cell.2014.06.037) 

![1-s2 0-S0092867414008642-gr4_lrg](https://github.com/user-attachments/assets/f5208b11-acc8-4473-9578-1bccd0c8e700)



## Calculating Beta Diversity Metrics

We will calculate the beta-diversity using the bray-curtis metric. 

First, let's normalize our dataset by relative abundance. 

```
phylo_obj_tree_sans_contam_sans_controls_relab <- transform_sample_counts(phylo_obj_tree_sans_contam_sans_controls, function(x) x / sum(x) )
```


```
bray_dist <- phyloseq::distance(phylo_obj_tree_sans_contam_sans_controls_relab, method="bray") # calculate bray-curtis metric
ordination <- ordinate(phylo_obj_tree_sans_contam_sans_controls_relab, method="PCoA", distance=bray_dist) # Perform ordination using bray-curtis metric
plot_ordination(phylo_obj_tree_sans_contam_sans_controls_relab, ordination, color="Site") + theme(aspect.ratio=1) # plot the ordination using ggplot2
```

We can get a list of the ordination methods available in phyloseq by requesting the help menu `??distanceMethodList`. Furthermore, we can test whether the sites differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:

```
adonis2(bray_dist ~ sample_data(phylo_obj_tree_sans_contam_sans_controls)$Site)
```
