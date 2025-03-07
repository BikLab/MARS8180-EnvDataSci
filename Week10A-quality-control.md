## The Metagenomic Dataset 
Starting today, we will start going to go through how to analyze metagenomic data. The test dataset, we are using here are 6 single-worm metagenomic sequences from two different families belonging to the family Oncholimidae and Epacanthion.

These worms were collected from Tybee Island in 2022 and 2023. They underwent REPLI-g amplification and sequencing using Illumina short-reads. Therefore, there are two major components in each metagenomic samples: 1) The host organisms and 2) the microbial members of the microbiome. 

We will go though the following steps through the next few modules: 

1. Quality Control of the Raw Metagenomic Data
2. Rapid Taxonomic Identification of Microbiome
2. Assembling Metagenomic Reads
3. Binning the contigs into Metagenome-assembled Genomes (MAGs)
4. Functional Annotations of metagenomic samples and MAGs 

Although we are analyzing a host-associated microbiome dataset, these same steps can be used to analyze any metagenomic dataset.

Before we begin, let's copy our files from the `instructors_data` directory. 

```
cp -r work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome /home/userid
```

You should have the following folders in your directory:

```
drwxr-xr-x. 2 ad14556 mars8180_instructor 16 Mar  3 16:56 data
drwxr-xr-x. 4 ad14556 mars8180_instructor  4 Mar  3 23:55 database
drwxr-xr-x. 2 ad14556 mars8180_instructor  2 Mar  3 15:06 metadata
drwxr-xr-x. 9 ad14556 mars8180_instructor 10 Mar  4 21:37 results
drwxr-xr-x. 2 ad14556 mars8180_instructor 10 Mar  4 21:59 scripts
```

## Visualizng Quality of Metagenomic Data

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
OUTPUT=//home/userid/nematode-microbiome/results/02-multiqc

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
