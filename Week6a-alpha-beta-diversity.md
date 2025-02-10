### Week 6a Alpha and Beta Diversity Analyses

#### Alpha Diversity

Some key points from Amy Willis's lectures at the STAMPS 2024 course at MBL: https://github.com/mblstamps/stamps2024/wiki#17

---

<img width="814" alt="Screenshot 2025-02-02 at 2 07 34 PM" src="https://github.com/user-attachments/assets/b737bb83-8771-4ee4-80e0-f1aade176d43" />
<img width="813" alt="Screenshot 2025-02-02 at 2 07 21 PM" src="https://github.com/user-attachments/assets/a8e990a8-3c84-4a68-9b46-973457236b5e" />

---

<img width="822" alt="Screenshot 2025-02-02 at 2 07 58 PM" src="https://github.com/user-attachments/assets/e7ee22c3-5c0c-405f-9de6-25b3421c7231" />

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
