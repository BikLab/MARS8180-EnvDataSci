![Untitled (6)](https://github.com/user-attachments/assets/e797dd79-0e91-4dcb-8d93-026604b0a5b7)
## The Metagenomic Dataset 

Starting today, we will start going to go through how to analyze metagenomic data. The test dataset we are using here are 6 single-worm metagenomic sequences belonging to two different families - Oncholaimidae and Thoracostomopsidae.

These worms were collected from Tybee Island in 2022 and 2023. They underwent REPLI-g amplification and sequencing using Illumina short-reads. Therefore, there are two major components in each metagenomic samples: 1) The host organisms and 2) the microbial members of the microbiome. 

We will go though the following steps through the next few modules: 

1. Quality Control of the Raw Metagenomic Data
2. Rapid Taxonomic Identification of the Microbiome
2. Assembling Metagenomic Reads
3. Binning the contigs into Metagenome-assembled Genomes (MAGs)
4. Functional Annotations of metagenomic samples and MAGs 

Although we are analyzing a host-associated microbiome dataset, these same steps can be used to analyze any metagenomic dataset.

Before we begin, let's copy our files from the `instructors_data` directory. 

```
cp -r /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome /home/userid
```

You should have the following folders in your directory:

```
drwxr-xr-x. 2 ad14556 mars8180_instructor 16 Mar  3 16:56 data
drwxr-xr-x. 4 ad14556 mars8180_instructor  4 Mar  3 23:55 database
drwxr-xr-x. 2 ad14556 mars8180_instructor  2 Mar  3 15:06 metadata
drwxr-xr-x. 9 ad14556 mars8180_instructor 10 Mar  4 21:37 results
drwxr-xr-x. 2 ad14556 mars8180_instructor 10 Mar  4 21:59 scripts
```

## Visualizing Quality of Metagenomic Data

As usual, we will start by going through the quality control of our metagenomic sequences. We will start by implementing `FASTQC` for every fastq file and collating that information with the `MULTIQC` package.

Let's first `cd` into our scripts subdirectory

```
cd /home/userid/nematode-microbiome/scripts/
```

Our first step is to run `FASTQC` and `MultiQC` on our raw metagenomic reads. 

Use nano to open and edit the first script `01-fastqc.sh`.

```
#!/bin/sh
#SBATCH --job-name="fastqc"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=user@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 01-fastqc.err-%N
#SBATCH -o 01-fastqc.out-%N

module load FastQC

INPUT=/home/userid/nematode-microbiome/data
OUTPUT=/home/userid/nematode-microbiome/results/01-fastqc

mkdir -p ${OUTPUT}
fastqc ${INPUT}/* -o ${OUTPUT} -t 8
```

Change the `userid` to your directory. This will run in parallel for every metagenomic sample in your `data` subdirectory. 

Afterwards let's edit the `multiqc` script. This is a tool that will collate all of our `fastqc` output file into one easily readable report. 

```
nano 02-multiqc.sh
```

```
#!/bin/sh
#SBATCH --job-name="multiqc"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 01-multiqc.err-%N
#SBATCH -o 01-multiqc.out-%N

module load MultiQC

INPUT=/home/userid/nematode-microbiome/results/01-fastqc
OUTPUT=/home/userid/nematode-microbiome/results/02-multiqc

multiqc --outdir ${OUTPUT} ${INPUT}
```

Submit these jobs to the cluster. 


## Running Trimmomatic
While these are running we can submit the trimmomatic script. First let's edit our script. 

```
nano 03-trimmomatic.sh
```

```
#!/bin/sh
#SBATCH --job-name="trimmomatic"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 03-trimmomatic.err-%N
#SBATCH -o 03-trimmomatics.out-%N

module load Trimmomatic

INPUT=/home/userid/nematode-microbiome/data
OUTPUT=/home/userid/nematode-microbiome/results/03-trimmomatic
ADAPTER=/home/userid/nematode-microbiome/database/trimmomatic/NexteraPE-PE.fa

mkdir ${OUTPUT}
for file in ${INPUT}/*_R1.fastq.gz; do
  base=$(basename ${file} _R1.fastq.gz)
  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 ${INPUT}/${base}_R1.fastq.gz ${INPUT}/${base}_R2.fastq.gz \
      ${OUTPUT}/${base}_R1_paired.fastq.gz ${OUTPUT}/${base}_R1_unpaired.fastq.gz \
      ${OUTPUT}/${base}_R2_paired.fastq.gz ${OUTPUT}/${base}_R2_unpaired.fastq.gz \
      SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${ADAPTER}:2:40:15 -threads 12
done
```

This script is more complicated but we will break it down. 

* First, we create a directory to save our results
* Second, we create a for loop to run this for each sample
* Afterwards, we use the `basename` function to get the sample names only
* Finally, we use the Trimmomatic module to trim our reads and remove the adapters. 

The flag `SLIDINGWINDOW:4:20` tells the program to trim whenever a group of 4 bases falls bellow a PHRED score of 20. `MINLEN:25` will only keep sequences with a minimum length of 25. The `ILLUMINACLIP:${ADAPTER}:2:40:15` will clip the adapters with 4 mismatches and a read score of 40 (highly recommend you read the manual). The adapter is located in this file `/home/userid/nematode-microbiome/database/trimmomatic/NexteraPE-PE.fa`

```
cat /home/userid/nematode-microbiome/database/trimmomatic/NexteraPE-PE.fa
```

```
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
```

We have the different combinations of the Nextera adapters, including the reverse complement sequences. 


## Visualize Quality of Trimmed Data

Now lets repeat the steps above to analyze the quality controlled sequences. Use nano to change `userid` to your 

```
nano 04-fastqc.sh 
nano 05-multqc.sh
```

## Viewing our MultiQC Report

To view our multiqc reports, first we need to download this to our computer. We are going to download our raw-sequence report and trimmed report to compare them.

```
scp -r ad14556@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/02-multiqc
```

```
scp -r ad14556@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/05-multiqc
```

The files will include an `html` file that we can open using any internet browser. 

## Taxonomic Profiling of Metagenomic Data

There are several tools you can use to assign taxonomy to short metagenomic reads, inlcuding [Kracken/Bracken](https://ccb.jhu.edu/software/bracken/), MetaPhlAn, [mOTUs2](https://www.nature.com/articles/s41467-019-08844-4), and [Centrifuge](https://ccb.jhu.edu/software/centrifuge/). Today, we will use MetaPhlAn4 to assign taxonomy to the nematode microbiome.

Description of MetaPhlAn from the Huttenhower Lab website (go here for user manual, tutorials, help forums, etc.): https://huttenhower.sph.harvard.edu/metaphlan/

MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. With StrainPhlAn, it is possible to perform accurate strain-level microbial profiling. MetaPhlAn 4 relies on ~5.1M unique clade-specific marker genes identified from ~1M microbial genomes (~236,600 references and 771,500 metagenomic assembled genomes) spanning 26,970 species-level genome bins (SGBs, http://segatalab.cibio.unitn.it/data/Pasolli_et_al.html), 4,992 of them taxonomically unidentified at the species level, allowing:

* unambiguous taxonomic assignments
* an accurate estimation of organismal relative abundance
* SGB-level resolution for bacteria, archaea and eukaryotes
* strain identification and tracking
* orders of magnitude speedups compared to existing methods.
* metagenomic strain-level population genomics

Software papers: 
 * **MetaPhlAn4**: Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., ... & Segata, N. (2023). [Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4](https://www.nature.com/articles/s41587-023-01688-w). _Nature Biotechnology_, 41(11), 1633-1644.
* **StrainPhlAn**: Truong, D. T., Tett, A., Pasolli, E., Huttenhower, C., & Segata, N. (2017). [Microbial strain-level population structure and genetic diversity from metagenomes](https://genome.cshlp.org/content/27/4/626.short). _Genome Research_, 27(4), 626-638.

<img width="926" alt="Screenshot 2025-03-13 at 8 51 00 AM" src="https://github.com/user-attachments/assets/74cff63b-a097-482f-9eb9-32d3a3122b05" />

(Figure 1 from Blanco-Míguez et al. 2023)

MetaPhlAn4 is more accurate than other software both in the detection of which taxa are present (higher F1 score) and their quantitative estimation (lower Bray-Curtis value asessing estimated taxon profiles vs. known abundance in gold standard dataset) - depicted Figure 2 from from Blanco-Míguez et al. 2023:

<img width="927" alt="Screenshot 2025-03-13 at 8 54 31 AM" src="https://github.com/user-attachments/assets/dceb783d-dee1-47bb-bb8e-a3f4727329dc" />

**However, note that every software paper will say their approach is the best! Authors assess their software using specific datasets (e.g. Human Microbiome Project studies) - You need to evaluate outputs from different tools on your own dataset!**

## Coding Tutorial

Let's `cd` into the scripts directory and use `nano` to edit the metaphlan4 script.


```
cd /home/userid/nematode-microbiome/scripts
nano 06A-metaphlan.sh
```

```
#!/bin/sh
#SBATCH --job-name="metaphlan"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 06-metaphlan.err-%N
#SBATCH -o 06-metaphlan.out-%N

module load Bowtie2
module load MetaPhlAn

INPUT=/home/userid/nematode-microbiome/results/03-trimmomatic
OUTPUT=/home/userid/nematode-microbiome/results/06-metaphlan
DATABASE=/home/userid/nematode-microbiome/database/metaphlan

mkdir -p ${OUTPUT}

for file in ${INPUT}/*_R1_paired.fastq.gz; do
  base=$(basename ${file} _R1_paired.fastq.gz)
  metaphlan ${INPUT}/${base}_R1_paired.fastq.gz,${INPUT}/${base}_R2_paired.fastq.gz --bowtie2out ${OUTPUT}/${base}_bowtie2.bz2 \
    --bowtie2db ${DATABASE} --ignore_eukaryotes --nproc 12 --input_type fastq -o ${OUTPUT}/${base}_taxonomy_profile.txt
done
```

Metphlan accepts paired-end sequences, but does not integrate paired-end data to assign taxonomy. We will use the trimmomatic trimmed paired-end sequences. We will also use the flag `--ignore_eukaryotes` to ignore any non bacterial sequences. We will run this in a loop for each sample and save them in the directory `/home/ad14556/nematode-microbiome/results/06-metaphlan`

After, we will collate that information using a helper script from Metaphlan4 and visualize them in a heatmap. 

```
nano 06B-metaphlan-viz.sh
```

```
#!/bin/sh
#SBATCH --job-name="metaphlan"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 06-metaphlan.err-%N
#SBATCH -o 06-metaphlan.out-%N

module load Bowtie2
module load MetaPhlAn

FOLDER=/home/userid/nematode-microbiome/results/06-metaphlan

merge_metaphlan_tables.py ${FOLDER}/*_profile.txt > ${FOLDER}/merged_abundance_table.txt
grep -E "g__|_taxonomy" ${FOLDER}/merged_abundance_table.txt | grep -v "s__" | sed "s/^.*|//g" | sed "s/_taxonomy//g" > ${FOLDER}/merged_abundance_table_genus.txt

hclust2.py \
-i ${FOLDER}/merged_abundance_table_genus.txt \
-o ${FOLDER}/merged_abundance_table_genus_heatmap.png \
--log_scale \
--dpi 30
```

First, we merge all of the metaphlan tables. After we use bash commands to data wrangle our data. 

* `grep -E "g__|_taxonomy"` get the lines that contain genus-level taxonomic information or samples
* `grep -v "s__"` removes lines that have species level information. 
* `sed "s/^.*|//g"` removes everything until the `g__` prefix
* `sed "s/_taxonomy//g"` removes `_taxonomy` from the sample names 

Finally, we will use the command hclust2.py to cluster and plot a heatmap of the log-transformed abundance table. 


We can download the heatmap to our computer to see the plot. 

```
scp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/06-metaphlan/merged_abundance_table_genus_heatmap.png
```

What bacteria might be members of the Oncholaimid core microbiome?

What bacteria might be members of the Epacanthion core microbiome?

## Taxonomic Profiling of Long-Read Metagenomic Data

Tools and algorithms for the taxonomoic profiling of metagenomic data were first developed for short-read sequences. These taxonomic profilers, such as metaphlan can NOT be used for long-read seqeunces because they utilize read alignment tools specifically designed for short-read sequences (<1,000 bp). Therefore, we are going to use a new tool, **Sourmash**, that is uses a k-mer approach to identify the taxonomic composition of unassembled long-read sequences. 

First, lets `cd` into our scripts directory and `nano` to edit our script

```
cd /home/ad14556/nematode-microbiome/scripts
nano 08-taxonomy-sourmash-long-reads.sh
```

```
#!/bin/bash

#SBATCH --job-name="sourmash"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e sourmash.err-%N
#SBATCH -o sourmash.out-%N

# path variables and modules
module load sourmash

INPUT=/home/userid/nematode-microbiome/data-long-read
OUTPUT=/home/userid/nematode-microbiome/results/08-taxonomy-long-read
DATABASE=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/database/sourmash

# single assembly
mkdir -p "$OUTPUT"

for FILE in "$INPUT"/*.fastq.gz; do
    SAMPLE=$(basename "$FILE" .fastq.gz)
    sourmash sketch dna -p scaled=1000,k=31 ${FILE} -o ${OUTPUT}/${SAMPLE}.sig
    sourmash gather ${OUTPUT}/${SAMPLE}.sig ${DATABASE}/gtdb-rs214-reps.k31.zip -o ${OUTPUT}/${SAMPLE}.csv --threshold-bp=3000
done

sourmash tax metagenome --gather-csv ${OUTPUT}/*.csv --taxonomy ${DATABASE}/gtdb-rs214.lineages.csv -o ${OUTPUT} --output-format kreport --rank genus
```

The sourmash developers have indexed various databases that can be directly used with their tool [https://sourmash.readthedocs.io/en/latest/databases.html](https://sourmash.readthedocs.io/en/latest/databases.html). To save time, we will be using the GTDB genomic representatives with 85,205 species-level genomes (indexed using a k-mer size of 31). For each each prepared database, they also made a file with the taxonomic information linking each genome with its assigned lineage (using either the GTDB or NCBI database as appropriate). **However, keep in mind that to accurately assign taxonomy to environmental samples you should use a more comprehensive database - in this case the Genbank bacterial database can be used to increase taxonomy assignments** This is especially important to K-mer approaches compared to alignment methods. K-mer methods look to **EXACT** sequence matches and a single SNP will cause that K-mer to remain unassigned.

The first step is to create a `sketch` or a list of k-mers that are seen in the metagenomic database. We are using a k-mer size of 31. A previous study has show that a k-mer size of either 31 or 51 have a high accuracy and precision in taxonomy assignments. The important part is the k-mer of the metagenome should match the database. 

Afterwards, the `gather` subcommand will be used to identify k-mers that are present in the reference genome. The `threshold-bp` indicates the the minimum estimated overlap for reporting a match. From the developers of sourmash:

> We have found a good intermediate threshold is 3 times the scaled value, e.g. --threshold-bp=3000 for a scaled value of 1000. This requires at least three overlapping hashes before a match is reported. If you are using a lower scaled value (a higher density sketch) because you are looking for matches between shorter sequences, then setting threshold-bp to 3 times that scaled value will take advantage of the increased sensitivity to short matches without introducing more false positives.

