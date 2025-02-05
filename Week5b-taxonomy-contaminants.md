### Week 5b Taxonomy Assignments & Alpha Diversity Analyses

#### Taxonomy Assignment


---
<img width="616" alt="Screenshot 2025-02-02 at 2 42 59 PM" src="https://github.com/user-attachments/assets/a9fae999-c335-47db-b64a-f5991fa657d4" />

Keck et al. 2022 (metabarcoding journal club paper): https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13746

---
From our lab (De Santiago et al, in revision at Environmental DNA):

Curation of taxonomy strings can be one of the most overlooked ways to improve reference databases:
<img width="1104" alt="Screenshot 2025-02-02 at 2 44 09 PM" src="https://github.com/user-attachments/assets/a94d6072-c136-4e7c-a7fc-c560806d4c11" />

ALEJANDRO TO INSERT CODE/EXERCISES HERE

---

#### Removing Contaminants

HOLLY TO INSERT THEORY / FIGURES / PAPER

### Assigning taxonomy to our sequences
Before we start, here are a list of important paths:

**Location of databases**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/database`

**Location of our scripts**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/scripts`

**Location of our results that I have run before class**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results
`

First, we need to copy over the script and database we are using to our personal home directory: 

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/database/*.qza \
  /home/userid/ddt-project/databases/
```

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/scripts/06-assign-taxonomy-blast.sh \
  /home/userid/ddt-project/scripts/
```

Let's cd into our scripts directory and nano into our `06-assign-taxonomy-blast.sh`

```
cd /home/userid/ddt-project/scripts/
nano 06-assign-taxonomy-blast.sh
```

```
#!/bin/sh
#SBATCH --job-name="assign-taxonomy"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=userid@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 06-assign-taxonomy.err-%N
#SBATCH -o 06-assign-taxonomy.out-%N

module load QIIME2/2024.10-amplicon

INPUT=/home/userid/ddt-project/results/05-dada2-rep-seq.qza
REFSEQ=/home/userid/ddt-project/databases/SILVA138-nr99_fixed_strings_custom_seq_sequences-May-16-2024.qza
REFTAX=/home/userid/ddt-project/databases/SILVA138-nr99_fixed_strings_custom_seq_taxonomy-May-16-2024.qza
CLASS=/home/userid/ddt-project/results/06-taxonomy-blast-90-1.qza
SEARCH=/home/userid/ddt-project/results/06-taxonomy-seach-results.qza

qiime feature-classifier classify-consensus-blast \
  --i-query ${INPUT} \
  --i-reference-taxonomy ${REFTAX} \
  --i-reference-reads ${REFSEQ} \
  --p-maxaccepts 1 \
  --p-perc-identity 0.90 \
  --o-classification ${CLASS} \
  --o-search-results ${SEARCH} \
  --p-num-threads 12
```

In this scripts, we are using the BLAST+ algorithm with the SILVA database with additional nematode sequences to assign taxonomic ID's to our ASVs. We will the top hit sequence with >90% sequence similarity. We will output the results in the file `06-taxonomy-blast-90-1.qza`.



First, we need to install and load required libraries to import and analyze.



We are going to read artifact files (.qza) using qiime2R


``` r
otus <- read_qza("/path/to/metadata/")
taxonomy <- read_qza("/path/to/taxonomy/")
tree <- read_qza("/path/to/tree/")
metadata <- read.delim2("/path/to/metdata/", sep = "\t", row.names = 1)
```

QZA files are zipped folders with many different pieces of information
including data provenance, format, version, and the data. We need to
extract out data from this file type

``` r
otu_df <- otus$data # we can view and save the actual data by specifying '$'
taxonomy_df <- taxonomy$data # save taxonomy info as a dataframe
phylo_tree <- tree$data # tree data is stored as a phyloseq object
```

The taxonomy file is not formatted correctly. All the taxonomy is in one
column. Ideally, each taxonomic level should have its own column.

``` r
head(taxonomy_df) 
```

    ##                         Feature.ID
    ## 1 00000ee244fdf5afdf2ffd124c23a320
    ## 2 0000a9d6f7dd5ee6cb6f2bed94784a39
    ## 3 00018b38663c5d40aa3c8e045f4d24aa
    ## 4 0002afd8dee5d314be3adc4cda9e46d2
    ## 5 0003175263e1b9530f1252de4f03e5ee
    ## 6 00033af273b8e729cf00555541979ad9
    ##                                                                                                                                                                        Taxon
    ## 1 D_0__Eukaryota;D_1__SAR;D_2__Alveolata;D_3__Ciliophora;D_4__Intramacronucleata;D_5__Conthreep;D_6__Plagiopylea;D_7__Odontostomatida;D_8__Epalxella;D_9__uncultured ciliate
    ## 2                                                D_0__Eukaryota;D_1__SAR;D_2__Rhizaria;D_3__Cercozoa;D_4__Vampyrellidae;D_5__uncultured;D_6__uncultured freshwater eukaryote
    ## 3  D_0__Eukaryota;D_1__SAR;D_2__Alveolata;D_3__Ciliophora;D_4__Intramacronucleata;D_5__Spirotrichea;D_6__Oligotrichia;D_7__Sinistrostrombidium;D_8__Strombidium paracalkinsi
    ## 4                                                                                                                                                                 Unassigned
    ## 5                                                            D_0__Eukaryota;D_1__SAR;D_2__Rhizaria;D_3__Cercozoa;D_4__Thecofilosea;D_5__uncultured;D_6__uncultured eukaryote
    ## 6                    D_0__Eukaryota;D_1__SAR;D_2__Stramenopiles;D_3__Ochrophyta;D_4__Diatomea;D_5__Bacillariophytina;D_6__Bacillariophyceae;D_7__uncultured marine eukaryote
    ##   Consensus
    ## 1         1
    ## 2         1
    ## 3         1
    ## 4         1
    ## 5         1
    ## 6         1

We are going to split the taxonomy into different columns and replace NA’s with with a placeholder: 'Unassigned'

``` r
taxonomy_fixed_df <- taxonomy_df %>% separate_wider_delim(Taxon, delim = ";", names_sep = "", too_few = "align_start")
taxonomy_fixed_df[is.na(taxonomy_fixed_df)] <- "Unassigned" # rename NAs into unassigned
head(taxonomy_fixed_df)
```

    ## # A tibble: 6 × 25
    ##   Feature.ID      Taxon1 Taxon2 Taxon3 Taxon4 Taxon5 Taxon6 Taxon7 Taxon8 Taxon9
    ##   <chr>           <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr> 
    ## 1 00000ee244fdf5… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… D_8__…
    ## 2 0000a9d6f7dd5e… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… Unass… Unass…
    ## 3 00018b38663c5d… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… D_8__…
    ## 4 0002afd8dee5d3… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unass…
    ## 5 0003175263e1b9… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… Unass… Unass…
    ## 6 00033af273b8e7… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… Unass…
    ## # ℹ 15 more variables: Taxon10 <chr>, Taxon11 <chr>, Taxon12 <chr>,
    ## #   Taxon13 <chr>, Taxon14 <chr>, Taxon15 <chr>, Taxon16 <chr>, Taxon17 <chr>,
    ## #   Taxon18 <chr>, Taxon19 <chr>, Taxon20 <chr>, Taxon21 <chr>, Taxon22 <chr>,
    ## #   Taxon23 <chr>, Consensus <dbl>

Now we can merge our otu table, taxonomy file, and tree into phylseq object

``` r
# fix taxonomy format 
taxonomy_fixed_df <- as.data.frame(taxonomy_fixed_df) # force into a dataframe
row.names(taxonomy_fixed_df) <- taxonomy_fixed_df$Feature.ID # make first column into row names
taxonomy_fixed_df$Feature.ID <- NULL # remove the first column
taxonomy_matrix <- as.matrix(taxonomy_fixed_df) # convert to a matrix 

physeq_otu <- otu_table(otu_df, taxa_are_rows = T) # convert into phyloseq object
physeq_tax <- tax_table(taxonomy_matrix) # convert into phyloseq object
physeq_meta <- sample_data(metadata) # convert into phyloseq object

phylo_object <- phyloseq(physeq_otu, physeq_tax, physeq_meta) # merge into phyloseq object
phylo_object_tree <- merge_phyloseq(phylo_object, phylo_tree) # add tree into phyloseq object
```

If we type in the phyloseq obect in our console, we can get a summary of our data. 

``` r
phylo_object_tree
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 190829 taxa and 4212 samples ]
    ## sample_data() Sample Data:       [ 4212 samples by 53 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 190829 taxa by 24 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 190829 tips and 190350 internal nodes ]



### Filter out contaminatns using Decontam’s prevelance method. For more information see <https://github.com/benjjneb/decontam>

As we've stressing throughout this course, we need to make sure we have a high-quality dataset. First, we can use our blanks to identify potential contaminations using a package called decontam. Let's set the prevalence theshold to to stricter value (0.5). For more information on paramters see the decontam manual

``` r
sample_data(phylo_object_tree)$is.neg <- sample_data(phylo_object_tree)$sample_control == "control" # create a sample-variable for contaminants
phylo_object_contaminants <- isContaminant(phylo_object_tree, method = "prevalence", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE) # detect contaminants based on control samples and their ASV prevalance
table(phylo_object_contaminants$contaminant) # check number of ASVs that are contaminents
```

    ## 
    ##  FALSE   TRUE 
    ## 190197    632

Now, we are going to make a presence-absence table of the contaminants in controls and samples. 

``` r
# Make phyloseq object of presence-absence in negative controls and true samples
phylo_object_contaminants.pa <- transform_sample_counts(phylo_object_tree, function(abund) 1 * (abund > 0)) # convert phyloseq table to presence-absence
ps.pa.neg <- prune_samples(sample_data(phylo_object_contaminants.pa)$sample_control == "control", phylo_object_contaminants.pa) # identify controls
ps.pa.pos <- prune_samples(sample_data(phylo_object_contaminants.pa)$sample_control == "true_sample", phylo_object_contaminants.pa) # identify samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=phylo_object_contaminants$contaminant) # convert into a dataframe
```

Let's plot our prevalance of ASVs in our blanks and compare then to our real samples. We see a clear split between prevalence in true samples vs controls.

``` r
# Make phyloseq object of presence-absence in negative controls and true samples
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```

![](01-filter-contamination-low-reads_files/figure-gfm/plot-contam-1.png)<!-- -->

Now, we can filter out these contaminants from our dataset.


``` r
phylo_obj_tree_sans_contam <- prune_taxa(!phylo_object_contaminants$contaminant, phylo_object_tree) # remove ASVs identified as decontaminants from the dataset
phylo_obj_tree_sans_contam_low <- filter_taxa(phylo_obj_tree_sans_contam, function(x) sum(x > 5) > 1, TRUE) # remove ASVs that are rare in each sample 
phylo_obj_tree_sans_contam_low_controls <- subset_samples(phylo_obj_tree_sans_contam_low, region != "kitblank" & region!= "negativecontrol" & region != "positivecontrol") ## Remove blanks and positive controls
phylo_obj_tree_sans_contam_low_controls
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 42455 taxa and 4095 samples ]
    ## sample_data() Sample Data:       [ 4095 samples by 54 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 42455 taxa by 24 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 42455 tips and 42385 internal nodes ]
