
##  Anvi'o

Anvi'o is a really powerful tool that allows you to analyze -omics datasets using modern bioinformatics tools. We are going to use it primarily to visualize the bins and manually assess and refine poorly binned contigs.

Before, we start there are a few terms specific to Anvi'o that you should know. 

###Contigs database

> A self-contained database containing a lot of information associated with your contig (or scaffold) sequences. This includes data that isn’t dependent on which sample the contigs came from, like positions of open reading frames, k-mer frequencies, split start/end points, functional and taxonomic annotations among others. You can initialize a basic contigs database from a FASTA file with the command anvi-gen-contigs-database and supplement it with additional information later in your analysis.

###Profile database

> A database containing sample-specific information about your contigs; for instance, coverage information from mapping reads to the contigs in a sample. Single profiles, each of which contains data for a particular sample, can be combined into a merged profile if they link to the same contigs database. The information across samples in a merged profile can be visualized as a ‘view’ in the anvi’o interactive database.

###Collection

> A virtual construct to store bins of items in an anvi’o profile database. Each collection contains one or more bins, and each bin contains one or more items. These items can be gene clusters, contigs, or other things depending on the display mode. See collection for more information.

## Installing Anvi'o

We will need to install Anvi'o using CONDA

```
interact --mem=15G

module load Miniconda3
mkdir /home/ad14556/conda-env/anvio
conda create -p /home/ad14556/conda-env/anvio
source activate /home/ad14556/conda-env/anvio

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

Since we simplified the contig names, we will have to rename them in the BAM files we previously made. Fortunately, the developers of anvi'o have a script for that.

```
anvi-script-reformat-bam sample.bam \
  -r contig-rename-report.txt \
  -o sample-rename.bam
```

Now we can create a sample profile for this database - this will include read mapping information. 

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

```
awk 'NR == FNR { a[$2] = $1; next } { $1 = a[$1] } 1' contig-rename-report.txt sample_DASTool_contig2bin.tsv > sample_DASTool_contig2bin-final.tsv
```
