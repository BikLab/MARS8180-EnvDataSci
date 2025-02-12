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

## Calculating Alpha Diversity
Now that we have created our phyloseq object, we are going to calculate Observed ASVS - a simple measure of alpha diversity. We can do this by using the `phyloseq::estimate_richness` function. We usually do not need to to specify the package, but this might be useful if you are using packages that have functions with the same name. 

```
alpha_div_observed <- phyloseq::estimate_richness(phylo_obj_tree_sans_contam_sans_controls, measures = "Observed") # calculate alpha diversity
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
metadata_df <- phyloseq::sample_data(phylo_obj_tree_sans_contam_sans_controls) # access metadata from the phyloseq object
head(metadata_df) 
```

We can put this together to calculate the observed diversity (number of ASVs) and create a dataframe that includes our metadata.

```
alpha_div_observed_metadata <- data.frame( # save as a dataframe
  phyloseq::sample_data(phylo_obj_tree_sans_contam_sans_controls), # access metadata 
  "Observed" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam_sans_controls, measures = "Observed")) # calculate Observed diversity and save it to a column names "Observed"
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
  phyloseq::sample_data(phylo_obj_tree_sans_contam_sans_controls), # get metadata
  "Reads" = phyloseq::sample_sums(phylo_obj_tree_sans_contam_sans_controls), # number of reads
  "Observed" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam_sans_controls, measures = "Observed"), # count observed ASVs
  "Shannon" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam_sans_controls, measures = "Shannon"), # Calculate Shannon Diversity
  "InvSimpson" = phyloseq::estimate_richness(phylo_obj_tree_sans_contam_sans_controls, measures = "InvSimpson")) # calculate InvSimpson
```

Now lets calculate Eveness by dividing Shannon by the Log of the Observed ASVs

```
alpha_div_observed_metadata$Evenness <- alpha_div_observed_metadata$Shannon/log(alpha_div_observed_metadata$Observed)
head(alpha_div_observed_metadata)
```

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
