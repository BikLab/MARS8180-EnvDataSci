## Read Mapping to our Assemblies

When we assemble our genome we lose quantitative information. We can map our reads back to the assembly to determine the relative abundance of each contig in the sample. This is important because we can use this information to bin our sequences into metagenome-assembled genomes (MAGs) - sequences that come from the same organisms should be present in roughly equal proportions. 

We will use two tools - **BWA** and **SAMTools** - to read map our samples. BWA will map short-reads to the assembly. After we will sort the BAM file and index it for rapid random access. 

* **BWA** (Li et al. 2009, Fast and accurate short read alignment with Burrows–Wheeler transform - [https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=false](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=false))
* **SamTools** (Danecek et al. 2021, Twelve years of SAMtools and BCFtools - [https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=false](https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=false))

To read map our samples, we 1) index our assembly, 2) read map using BWA, 3) sort to SAM file and convert to BAM, and 4) index the sorted BAM file.

In practice it looks like this:

```
bwa index contig_file.fa 
bwa mem contig_file.fa forward_reads.fastq reverse_reads.fastq -o alignment.sam
samtools sort alignment.sam -o alignment.bam
samtools index alignment.bam
``` 

Lets first copy the scripts to our directory

```
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/10-read-mapping-short-reads.sh /home/userid/nematode-microbiome/scripts
```

Now use the text editor to edit the script: 

```
#!/bin/bash

#SBATCH --job-name="read-map"
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e read-map.err-%N
#SBATCH -o read-map.out-%N

# path variables and modules
module load BWA
module load SAMtools

CONTIGS=/scratch/userid/nematode-microbiome/results/07-assembly
READS=/home/userid/nematode-microbiome/results/03-trimmomatic
OUTPUT=/scratch/userid/nematode-microbiome/results/10-read-map-short-reads

for FILE in ${CONTIGS}/*; do
  SAMPLE=$(basename ${FILE})
  mkdir -p ${OUTPUT}/${SAMPLE}/
  bwa index ${CONTIGS}/${SAMPLE}/final.contigs.fa
  bwa mem -t 48 ${CONTIGS}/${SAMPLE}/final.contigs.fa ${READS}/${SAMPLE}_R1_paired.fastq.gz ${READS}/${SAMPLE}_R2_paired.fastq.gz | samtools sort -o ${OUTPUT}/${SAMPLE}/${SAMPLE}-alignment.bam --threads 48
  samtools index -@ 48 ${OUTPUT}/${SAMPLE}/${SAMPLE}-alignment.bam
done
```

## Metagenome-assembled Genomes

We are now able to use our abundance information to bin bacterial contigs into metagenome-assembled genomes. We are going to use three tools 1) MetaBat2, 2) comebin, and 3) dastool. Both MetaBat2 and Comebin use tetranucleotide frequencies in conjunction with abundance information for genome reconstruction. However, Comebin uses a contrastive multi-view representation learning to determine the best MAGs and incorporates single-copy gene information and contig length. Finally, Dastool compares the MAGs produced by any binning algorithm and determine the most complete MAGs (with least amount of contamination). 

* **MetaBat2** (Kang et al. 2019, MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies - [https://peerj.com/articles/7359/](https://peerj.com/articles/7359/))
* **Comebin** (Wang et al. 2024, Effective binning of metagenomic contigs using contrastive multi-view representation learning - [https://www.nature.com/articles/s41467-023-44290-z](https://www.nature.com/articles/s41467-023-44290-z))
* **Dastool** (Sieber et al. 2018, Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy - [https://www.nature.com/articles/s41564-018-0171-1](https://www.nature.com/articles/s41564-018-0171-1))

### MetaBat2

![peerj-03-1165-g001](https://github.com/user-attachments/assets/d3604d96-0238-42dd-a283-4590f7fae801)

The figure above is from the MetaBat1 software; howwever, MetaBat2 works similarly, except that it is iterative and chooses the best parameters according to your dataset. It does not require to make any decisions about which parameters would be best suited for your dataset. Additionally, MetaBat2 constructs a graph using the tetranucleotide frequency and read abundance to cluster contigs that appear to be similar in structure. 

For MetaBat2, we are first going to summarize our BAM file and then reconstruct the genomes.

```
jgi_summarize_bam_contig_depths --outputDepth contig-depth.txt alignment.bam
metabat2 -i final.contigs.fa -a contig_depth -o sample-name
```

### Comebin

![41467_2023_44290_Fig6_HTML](https://github.com/user-attachments/assets/fa8c5aee-699a-47e3-aeb1-0e91badc0ded)

From the manuscript, "The key contributions of COMEBin can be summarized as follows: 1) We introduce a data augmentation approach that generates multiple views for each contig, enabling contrastive learning and yielding high-quality representations of the heterogeneous features; 2) COMEBin incorporates a “Coverage module” to obtain fixed-dimensional coverage embeddings, which enhances its performance across datasets with varying numbers of sequencing samples; 3) COMEBin employs the advanced community detection algorithm, Leiden21, for clustering. Moreover, we adapt the settings of Leiden specifically for the binning task, considering single-copy gene information22 and contig length. This adaptation ensures that COMEBin produces robust and reliable binning results across diverse datasets." 

![41467_2023_44290_Fig2_HTML](https://github.com/user-attachments/assets/4f6435c3-a8c4-4d25-b060-d35cedf98fd4)

Comebin has been shown to outperform every popular binning algorithm and results in MAGs with higher genome completeness and least amount of contamination (more on this on Thursday). 

Comebin can be run with a single-line

```
run_comebin.sh -a final.contigs.fa -o sample-directory -p alignment.bam -t 40
```

### DAStool 

![image](https://github.com/user-attachments/assets/d6cd0bbb-ed5c-42df-8625-76b8963998d0)

* Step 1: The DasTool inputs are the bins output by each binning algorith. 
* Step 2: Single-copy genes are predicted and the bins are scored according to completeness and contamination.
* Step 3: Bins that were predicted by multiple binning algorithms are dereplicated.
* Step 4: There is an iterative selection of "best" bins and the partial candidate bins are selected. The output are a non-redundant set of the high-scoring bins from predicted by different binning algorithms.

To run DASTool we first need to make a list of contigs that belong to each bin. Afterwards we can compare the assemblies and choose the best one. We can set the score-threshold to 0 to force it to bin incomplete MAGs (low completion according to single-copy genes). 

```
Fasta_to_Contig2Bin.sh -i metabat-bins -e fa > metabat-summary.txt
Fasta_to_Contig2Bin.sh -i comebin-bins -e fa > comebin-summary.txt
Rscript DAS_Tool.R -i comebin-summary.txt,metabat-summary.txt -l comebin,metabat -c final.contigs.fa -o dastool-bins --write_bins --write_bin_evals -t 12 --score_threshold=0
```

### In Practice: Binning our Nematode-associated Metagenome-Assembled Genomes

Let's copy the binning scripts to our home directory.

```
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/11-metabat2-short-reads.sh /home/userid/nematode-microbiome/scripts 
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/12-comebin-short-reads.sh /home/userid/nematode-microbiome/scripts 
```

Let's edit the metabat2 using a text editor. 

```
cd /home/userid/nematode-microbiome/scripts 
nano 11-metabat2-short-reads.sh
```

```
#!/bin/bash

#SBATCH --job-name="metabat"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=userid@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e metabat.err-%N
#SBATCH -o metabat.out-%N

# path variables and modules
module load MetaBAT

CONTIGS=/scratch/ad14556/nematode-microbiome/results/07-assembly
MAP=/scratch/userid/nematode-microbiome/results/10-read-map-short-reads
DEPTH=/scratch/userid/nematode-microbiome/results/11-metabat-short-reads

for FILE in ${CONTIGS}/*; do
  mkdir ${DEPTH}/${SAMPLE}
  SAMPLE=$(basename ${FILE})
  jgi_summarize_bam_contig_depths --outputDepth ${DEPTH}/${SAMPLE}-depth.txt ${MAP}/${SAMPLE}/${SAMPLE}-alignment.bam
  metabat2 -i ${CONTIGS}/${SAMPLE}/final.contigs.fa -a ${DEPTH}/${SAMPLE}-depth.txt -o ${DEPTH}/${SAMPLE}/${SAMPLE} -t 24
done
```
For comebin, we will first need to install this as a conda environment since it is not currently installed on the cluster

```
interact --partition=batch --mem=15G
module load Miniconda3

mkdir /home/ad14556/conda-env/comebin
conda create -p /home/ad14556/conda-env/comebin
source activate /home/ad14556/conda-env/comebin
conda install -c conda-forge -c bioconda comebin
```

While this is installing, we can create/edit our script to use this software to bin our metagenomes. 

```
nano 12-comebin-short-reads.sh
```

```
#!/bin/bash

#SBATCH --job-name="comebin"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=userid@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e comebin.err-%N
#SBATCH -o comebin.out-%N

# path variables and modules
module load MetaBAT

CONTIGS=/scratch/ad14556/nematode-microbiome/results/07-assembly
MAP=/scratch/userid/nematode-microbiome/results/10-read-map-short-reads
BIN=/scratch/userid/nematode-microbiome/results/12-comebin-short-reads

for FILE in ${CONTIGS}/*; do
  mkdir ${BIN}/${SAMPLE}
  SAMPLE=$(basename ${FILE})
  run_comebin.sh -a ${CONTIGS}/${SAMPLE}/final.contigs.fa -o ${BIN}/${SAMPLE} -p ${MAP}/${SAMPLE} -t 24
done
```
