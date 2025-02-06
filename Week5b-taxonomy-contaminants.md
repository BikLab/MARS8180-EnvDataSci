### Week 5b Taxonomy Assignments & Alpha Diversity Analyses

#### Taxonomy Assignment


---
<img width="616" alt="Screenshot 2025-02-02 at 2 42 59 PM" src="https://github.com/user-attachments/assets/a9fae999-c335-47db-b64a-f5991fa657d4" />

Keck et al. 2022 (metabarcoding journal club paper): https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13746

---
From our lab (De Santiago et al, in revision at Environmental DNA):

Curation of taxonomy strings can be one of the most overlooked ways to improve reference databases:
<img width="1104" alt="Screenshot 2025-02-02 at 2 44 09 PM" src="https://github.com/user-attachments/assets/a94d6072-c136-4e7c-a7fc-c560806d4c11" />

---

Hleap et al. (2021) Assessment of current taxonomic assignment strategies for metabarcoding eukaryotes, Molecular Ecology Resources - https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13407

Cortina et al (2023) Improving species-level taxonomic assignments from 16S rRNA sequencing technologies - https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.930 

---

Four main taxonomy assignment methods:
* Sequence similarity (SS)
* Sequence composition (SC)
* Phylogenetic methods (Ph)
* Probabalistic methods (Pr)

<img width="1137" alt="Screenshot 2025-02-06 at 9 06 16 AM" src="https://github.com/user-attachments/assets/52636e4f-4629-4941-9382-ed6add417797" />


#### Removing Contaminants

HOLLY TO INSERT THEORY / FIGURES / PAPER

## Assigning taxonomy to our sequences
Before we start, here are a list of important paths:

**Location of metadata**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/metadata`

**Location of databases**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/database`

**Location of our scripts**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/scripts`

**Location of our results**: `/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results
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


## Downloading our data to our personal computer

Before we download our ASV table, taxonomy table, and tree onto our personal computer, we are going to create a projects directory on our personal laptop. On your desktop, create a folder called `ddt-project` where we will store our QIIME2 files and our downstream data analysis. 

The file path will vary for each user and computing system, so make sure you are using the file path on your computer. 

```
cd /Users/userid/Desktop
mkdir ddt-project
cd ddt-project
mkdir results data scripts metadata
```

Now, we can download our data using the scp command.

**Download the metadata (`2025-01-03-ddt-metadata.csv`)** 

```
scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/metadata/2025-01-03-ddt-metadata.csv metadata/
```

**Download the ASV Table (`05-dada2-feature-table.qza`)** 

```
scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results/05-dada2-feature-table.qza results/
```

**Download the Taxonomy Table (`06-taxonomy-blast-90-1.qza`)** 

```
scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results/06-taxonomy-blast-90-1.qza results/
```

**Download the Phylogenetic Tree (`07-fasttree-midrooted-tree.qza`)** 

```
scp userid@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results/07-fasttree-midrooted-tree.qza results/
```

## Starting an R project

Now, we can start an R project in the directory we just created.

Let's install and import the R packages we are going to use for our downstream analysis. **QIIME2R** allows us to easily import the artifact files into a "phyloseq" object. **Phyloseq** let's us manipulate our metabarcoding dataset. **Decontam** is a package that allows us to use our blank to remove potential contaminants. **Tidyr** and **ggplot** allow us to easily manipulate dataframes and create publication ready plots, respectively.

```
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}

devtools::install_github("jbisanz/qiime2R") # install qiime2R 
BiocManager::install("decontam") # install decontam
BiocManager::install("phyloseq") # install phyloseq
BiocManager::install("ggplot2") # install ggplot2

library(phyloseq)
library(decontam)
library(qiime2R)
library(tidyr)
library(ggplot2)
```

## Importing our data into phyloseq object
We are going to read artifact files (.qza) using qiime2R and import our metadata file using the `read.delim2` command.


``` r
otus <- read_qza("results/05-dada2-feature-table.qza")
taxonomy <- read_qza("results/06-taxonomy-blast-90-1.qza")
tree <- read_qza("results/07-fasttree-midrooted-tree.qza")
metadata <- read.delim2("metadata/2025-01-03-ddt-metadata.csv", sep = ",", row.names = 1)
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

```
                        Feature.ID
1 0005b35778658ed77217967b57fdd319
2 0006630fb18ffea20c1bd1c0227f22f3
3 000fb3ff10896c00c9c11d750da7a164
4 0013dd4bef6114837fdf8c29ee55ea50
5 0017c06698f0b928e93ffd3fd8498011
6 001b7c34dad08db30819fdf95a704b3e
                                                                                                                                                                                                                    Taxon
1                                                                                                      D_0__Eukaryota;D_1__Amorphea;D_2__Amoebozoa;D_3__Incertae Sedis;D_4__Apusomonadidae;D_5__uncultured Apusomonadidae
2                                                                                                                D_0__Eukaryota;D_1__SAR;D_2__Rhizaria;D_3__Cercozoa;D_4__Novel Clade 12;D_5__uncultured marine eukaryote
3                                                                                                 D_0__Eukaryota;D_1__SAR;D_2__Rhizaria;D_3__Retaria;D_4__Polycystinea;D_5__Collodaria;D_6__AT8-54;D_7__marine metagenome
4 D_0__Eukaryota;D_1__Amorphea;D_2__Obazoa;D_3__Opisthokonta;D_4__Nucletmycea;D_5__Fungi;D_6__Dikarya;D_7__Basidiomycota;D_8__Agaricomycotina;D_9__Agaricomycetes;D_10__Agaricales;D_11__Chamaeota;D_12__Chamaeota sinica
5                                                                                           D_0__Eukaryota;D_1__SAR;D_2__Alveolata;D_3__Protalveolata;D_4__Syndiniales;D_5__Syndiniales Group I;D_6__uncultured eukaryote
6                                                                                                                                                                                                              Unassigned
  Consensus
1         1
2         1
3         1
4         1
5         1
6         1
```

We are going to split the taxonomy into different columns and replace NA’s with with a placeholder: 'Unassigned'

``` r
taxonomy_fixed_df <- taxonomy_df %>% separate_wider_delim(Taxon, delim = ";", names_sep = "", too_few = "align_start")
taxonomy_fixed_df[is.na(taxonomy_fixed_df)] <- "Unassigned" # rename NAs into unassigned
head(taxonomy_fixed_df)
```

```
# A tibble: 6 × 25
  Feature.ID       Taxon1 Taxon2 Taxon3 Taxon4 Taxon5 Taxon6 Taxon7 Taxon8 Taxon9 Taxon10 Taxon11
  <chr>            <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>   <chr>  
1 0005b35778658ed… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… Unass… Unass… Unass… Unassi… Unassi…
2 0006630fb18ffea… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… Unass… Unass… Unass… Unassi… Unassi…
3 000fb3ff10896c0… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… Unass… Unassi… Unassi…
4 0013dd4bef61148… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… D_8__… D_9__A… D_10__…
5 0017c06698f0b92… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… Unass… Unass… Unassi… Unassi…
6 001b7c34dad08db… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unassi… Unassi…
# ℹ 13 more variables: Taxon12 <chr>, Taxon13 <chr>, Taxon14 <chr>, Taxon15 <chr>,
#   Taxon16 <chr>, Taxon17 <chr>, Taxon18 <chr>, Taxon19 <chr>, Taxon20 <chr>, Taxon21 <chr>,
#   Taxon22 <chr>, Taxon23 <chr>, Consensus <dbl>
```

Let's force our taxonomy table back into a matrix datatype. 

``` r
taxonomy_fixed_df <- as.data.frame(taxonomy_fixed_df) # force into a dataframe
row.names(taxonomy_fixed_df) <- taxonomy_fixed_df$Feature.ID # make first column into row names
taxonomy_fixed_df$Feature.ID <- NULL # remove the first column
taxonomy_matrix <- as.matrix(taxonomy_fixed_df) # convert to a matrix 
```

Now we can merge our otu table, taxonomy file, and tree into a phyloseq object
```
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

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24188 taxa and 208 samples ]
sample_data() Sample Data:       [ 208 samples by 43 sample variables ]
tax_table()   Taxonomy Table:    [ 24188 taxa by 24 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 24188 tips and 24167 internal nodes ]
```
---
#### Removing Contaminant Sequences

<img width="762" alt="Screenshot 2025-02-06 at 9 09 16 AM" src="https://github.com/user-attachments/assets/71a216ca-2e93-4a71-aa68-a555bcf2e50a" />

Image from Davis et al. (2018) Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data, Microbiome - https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2

## Filter out contaminatns using Decontam’s prevelance method. For more information see <https://github.com/benjjneb/decontam>

As we've stressing throughout this course, we need to make sure we have a high-quality dataset. First, we can use our blanks to identify potential contaminations using a package called decontam. Let's set the prevalence theshold to to stricter value (0.5). For more information on parameters see the decontam manual

``` r
sample_data(phylo_object_tree)$is.neg <- sample_data(phylo_object_tree)$Sample_Control == "Control" # create a sample-variable for contaminants
phylo_object_contaminants <- isContaminant(phylo_object_tree, method = "prevalence", neg="is.neg", threshold=0.6, detailed = TRUE, normalize = TRUE) # detect contaminants based on control samples and their ASV prevalance
table(phylo_object_contaminants$contaminant) # check number of ASVs that are contaminents
```
``` 
##  FALSE   TRUE 
## 190197    632
```
Now, we are going to make a presence-absence table of the contaminants in controls and samples. 

``` r
# Make phyloseq object of presence-absence in negative controls and true samples
phylo_object_contaminants.pa <- transform_sample_counts(phylo_object_tree, function(abund) 1 * (abund > 0)) # convert phyloseq table to presence-absence
ps.pa.neg <- subset_samples(phylo_object_contaminants.pa, Sample_Control=="Control")
ps.pa.pos <- subset_samples(phylo_object_contaminants.pa, Sample_Control=="Sample")
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=phylo_object_contaminants$contaminant) # convert into a dataframe
```

Let's plot our prevalance of ASVs in our blanks and compare then to our real samples. We see a clear split between prevalence in true samples vs controls.

``` r
# Make phyloseq object of presence-absence in negative controls and true samples
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```

Now, we can filter out these contaminants from our dataset.


``` r
phylo_obj_tree_sans_contam <- prune_taxa(!phylo_object_contaminants$contaminant, phylo_object_tree) # remove ASVs identified as decontaminants from the dataset
phylo_obj_tree_sans_contam_sans_controls <- subset_samples(phylo_obj_tree_sans_contam, Sample_Control != "Control") ## Remove blanks and positive controls
phylo_obj_tree_sans_contam_sans_controls
```

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 16568 taxa and 202 samples ]
sample_data() Sample Data:       [ 202 samples by 43 sample variables ]
tax_table()   Taxonomy Table:    [ 16568 taxa by 24 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 16568 tips and 16563 internal nodes ]
```
