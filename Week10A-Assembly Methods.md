## Assembling Short-Read Sequences

There are several tools you can use to assembly short-read metagenomics, but we are going to implement MegaHit for three main reasons: 1) it is incredibly fast, 2) it requires less computational resources, and 3) assembles metagenomic datasets fairly well. 

First, lets `cd` into the scripts directory: 

```
/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts
```

Now, let's use `nano` to edit the assembly script. 

```
nano 07-megahit.sh
```

```
#!/bin/sh
#SBATCH --job-name="megahit"
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=30G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 07-megahit.err-%N
#SBATCH -o 07-megahit.out-%N

module load MEGAHIT

INPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/03-trimmomatic
OUTPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/07-assembly

mkdir ${OUTPUT}

for file in ${INPUT}/*_R1_paired.fastq.gz; do 
  base=$(basename ${file} _R1_paired.fastq.gz)
  megahit -1 ${INPUT}/${base}_R1_paired.fastq.gz -2 ${INPUT}/${base}_R2_paired.fastq.gz -o ${OUTPUT}/${base} -t 24 --presets meta-sensitive
done
```

We are going to assemble each sample individually. Although most of these nematodes were collected from Tybee Island, they were isolated from different habitats which can affect the variability and composition of the microbiome. 

We have to first create a for loop for each sample. Afterwards we save the `basename` of each sample. Afterwards we can use the `megahit` software to assemble our paired-end reads (R1 and R2). We will use a preset called `meta-sensitive` which specifies the list of Kmers to use to build succint de brujin graphs. 

**Important Definition**

* **Kmers:** A substring of a sequence of the length K 
	
	Take the sequence ATGCTGCT as an example. This sequence can be broken up into several substrings (or subsequences) of length 4 (tetramer or 4-mer).
	* ATGC
	* TGCT
	* GCTG
	* CTGC
	* TGCT

	This is useful because we can count how many unique and total kmers we find in a sequence and quickly compare them to other sequences. In the example above, there are 5 tetramers (4 unique tetramers). This is useful because alignment algorithms require a lot of resources - the time and space required grow exponentially for each sequence. Using the same example above, there are only 256 (4^4) unique tetramers. If we increase the size of the kmer, we have more unique kmer combination but fewer overlapping substrings. Even BLAST utilizes Kmers (default size 28; 4^28 unique combinations) to identify exact matches across the entire NCBI database (>3.7 billion sequences). 

* **De brujin Graphs:** Graph theory underpin many -omics assembly methods. De brujin graphs are old - they were first developed in 1946 by the Mathmatecian Nicolaas de Brujin. In short they are a directed graph-based method of visualizing and assembling sequence data. You take your kmers and connect them if they overlap. Afterwards, you can follow all your overlapping sequences to identify assembled contigs. 

The megahit `meta-sensitive` uses a kmer list of `21,31,41,51,61,71,81,91,99` to create debrujin graphs. Afterwards, it merges all the graphs together. 
