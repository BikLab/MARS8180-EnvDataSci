
##  Anvi'o

Anvi'o is a really powerful tool that allows you to analyze -omics datasets using modern bioinformatics tools. We are going to use it primarily to visualize the bins and manually assess and refine poorly binned contigs.

Before, we start there are a few terms specific to Anvi'o that you should know. 

### Contigs database

> A self-contained database containing a lot of information associated with your contig (or scaffold) sequences. This includes data that isn’t dependent on which sample the contigs came from, like positions of open reading frames, k-mer frequencies, split start/end points, functional and taxonomic annotations among others. You can initialize a basic contigs database from a FASTA file with the command anvi-gen-contigs-database and supplement it with additional information later in your analysis.

### Profile database

> A database containing sample-specific information about your contigs; for instance, coverage information from mapping reads to the contigs in a sample. Single profiles, each of which contains data for a particular sample, can be combined into a merged profile if they link to the same contigs database. The information across samples in a merged profile can be visualized as a ‘view’ in the anvi’o interactive database.

### Collection

> A virtual construct to store bins of items in an anvi’o profile database. Each collection contains one or more bins, and each bin contains one or more items. These items can be gene clusters, contigs, or other things depending on the display mode. See collection for more information.

## Installing Anvi'o

We will need to install Anvi'o using CONDA

```
interact --mem=15G

module load Miniconda3
mkdir /home/userid/conda-env/anvio
conda create -p /home/userid/conda-env/anvio
source activate /home/userid/conda-env/anvio

conda install -y -c conda-forge -c bioconda python=3.10 \
        sqlite=3.46 prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript \
        nodejs=20.12.2
conda install -y -c bioconda pysam

curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
        --output anvio-8.tar.gz
pip install anvio-8.tar.gz
```

There are several steps to do this. First, we need to reformat the sequence names in our fasta file. We can do this using the command `anvi-script-reformat-fasta fasta`. This will create a file with simplified contig names and save the changes in a report file. We will need this later. 

```
anvi-script-reformat-fasta contigs.fa \
	-o contigs-renamed.fa \
	--simplify-names \
	--report-file contig-rename-report.txt
```

Afterwards, we can create and anvi'o database using `anvi-gen-contigs-database`. This command can call on several tools to annotate your data including: 

* prodigal (Gene calling)
* HMMER (HMM search)
* krakenuniq (Gene taxonomy)
* centrifuge (Gene taxonomy)
* DIAMOND (Sequence search against various databases)

```
anvi-gen-contigs-database -f contigs-renamed.fa \
	-o sample.db -n sample-name 
```

Since we simplified the contig names, we have to make sure the contig names from the BAM files match. **We will have to redo the read mapping**. Afterwards, we can create a sample profile for this database - this will include read mapping information. 

```
anvi-profile -i  sample-rename.bam \
	-c sample.db \
	-o profile.db
```

Finally, we want to add our bins into a collection. But if we look at our contig-bin mapping information, we still have old names associated with our contigs. 

```
head /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/13-dastool-short-reads/epacanthion.1/epacanthion_003_GA_TI_202311_DASTool_contig2bin.tsv
```

```
k141_1003579	10630
k141_1010370	10630
k141_1012047	10630
k141_1028068	10630
k141_1038723	10630
k141_103904	10630
k141_1039738	10630
k141_1044728	10630
k141_1046824	10630
k141_1056006	10630
```

So, we are going to have to change that using a bash commands - we will need the report file we created earlier. 

**I have already run the commands below for sample Epacanthion.1**

```
# change contig names to match changes made by anvio
awk 'NR == FNR { a[$2] = $1; next } { $1 = a[$1] } 1' contig-rename-report.txt sample_DASTool_contig2bin.tsv > sample_DASTool_contig2bin-final.tsv

# replace spaces with tabs
sed 's/  \+/\t/g' sample_DASTool_contig2bin-final.tsv > dastool-scaffolds2bin-tab-final-no-spaces.tsv

# add sample name to bins (anvio doesnt like it when names start with numebers)
sed 's/\t/\t\epacanthion-1-/' dastool-scaffolds2bin-tab-final-no-spaces.tsv >> dastool-scaffolds2bin-tab-final-no-spaces-rename-bins.tsv

```

Afterwards, we can add our bins to the collection: 

```
anvi-import-collection sample-dastool-scaffolds2bin-tab-final-no-spaces-rename-bins.tsv -C description -p PROFILE -c sample.db --contigs-mode
```

Let's copy the scripts from the instructors directory.

```
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/17A-anvio-read-mapping.sh /home/userid/nematode-microbiome/scripts
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/17B-anvio-viz.sh /home/userid/nematode-microbiome/scripts
```

First, we need to read map one more time so that the file has the correct contigs names

```
nano /home/userid/nematode-microbiome/scripts/17A-anvio-read-mapping.sh
```

```
#!/bin/bash

#SBATCH --job-name="read-map"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e read-map.err-%N
#SBATCH -o read-map.out-%N

# path variables and modules
module load BWA
module load SAMtools

CONTIGS=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/18-anvio-short-reads/epacanthion.1/final-contigs-reformat.fa
FORWARD=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/03-trimmomatic/epacanthion.1_R1_paired.fastq.gz
REVERSE=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/03-trimmomatic/epacanthion.1_R2_paired.fastq.gz
OUTPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/18-anvio-short-reads/epacanthion.1

bwa index ${CONTIGS}
bwa mem -t 24 ${CONTIGS} ${FORWARD} ${REVERSE} | samtools sort -o ${OUTPUT}/epacanthion.1-alignment-rename.bam --threads 48
samtools index -@ 24 ${OUTPUT}/epacanthion.1-alignment-rename.bam
```

Then, we can run through the Anvi'o workflow:

```
nano /home/userid/nematode-microbiome/scripts/17B-anvio-viz.sh
```

```
#!/bin/bash

#SBATCH --job-name="anvio"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e anvio.err-%N
#SBATCH -o anvio.out-%N

# path variables and modules
module load Miniconda3
source activate /home/ad14556/conda-env/anvio

CONTIGS=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/07-assembly/epacanthion.1/final.contigs.fa
BINS=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/13-dastool-short-reads/epacanthion.1/epacanthion_003_GA_TI_202311_DASTool_contig2bin.tsv
BAM=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/18-anvio-short-reads/epacanthion.1/epacanthion.1-alignment-rename.bam
OUTPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/18-anvio-short-reads/epacanthion.1

# simplify contig names
anvi-script-reformat-fasta ${CONTIGS} -o ${OUTPUT}/final-contigs-reformat.fa --simplify-names --report-file ${OUTPUT}/final-contigs-reformat.txt
# create the contigs database
anvi-gen-contigs-database -f ${OUTPUT}/final-contigs-reformat.fa -o ${OUTPUT}/epacanthion.1.db -n epacanthion.1 -T 24
# create single sample profile
anvi-profile -i ${OUTPUT} -c ${OUTPUT}/epacanthion.1.db -o ${OUTPUT}/profile --sample-name epacanthion_1 -T 24
#import bins
anvi-import-collection ${OUTPUT}/epacanthion-1-dastool-scaffolds2bin-tab-final-no-spaces-rename-bins.tsv -C dastool -p ${OUTPUT}/profile/PROFILE.db -c ${OUTPUT}/epacanthion.1.db --contigs-mode
```

## METABOLIC

There are other tools that are not available via Anvi'o that allow you to predict metabolic and biogeochemical functional trait from metagenomic data. One of the tools is METABOLIC. The genomic datasets can either be metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or isolated strain sequenced genomes. There are two programs we can run with METABOLIC - METABOLIC-C and METABOLIC-G. 

GitHub: [https://github.com/AnantharamanLab/METABOLIC](https://github.com/AnantharamanLab/METABOLIC)

### METABOLIC-G.pl	
> Allows for classification of the metabolic capabilities of input genomes.

### METABOLIC-C.pl	
> Allows for classification of the metabolic capabilities of input genomes,
> calculation of genome coverage, creation of biogeochemical cycling diagrams,
> and visualization of community metabolic interactions and contribution to biogeochemical processes by each microbial group.

For our purposed, we are going to use METABOLIC-G.  

```
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/18-metabolic.sh /home/userid/nematode-microbiome/scripts/17A-anvio-read-mapping.sh
```

```
nano /home/userid/nematode-microbiome/scripts/17A-anvio-read-mapping.sh
```

```
#!/bin/sh

#SBATCH --job-name="metabolic"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e metabolic.err-%N
#SBATCH -o metabolic.out-%N

#Path Variables
module load HMMER
module load METABOLIC
module load GTDB-Tk

INPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/16-final-bins-short-reads/epacanthion.1
OUTPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/19-metabolic/epacanthion.1

mkdir -p ${OUTPUT}/input
mkdir -p ${OUTPUT}/output

cp ${INPUT}/*.fa ${OUTPUT}/input
rename ".fa" ".fasta" ${OUTPUT}/input/*

METABOLIC-G.pl -p meta -t 24 -in-gn ${OUTPUT}/input -o ${OUTPUT}/output
```
