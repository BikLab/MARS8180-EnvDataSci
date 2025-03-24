## Read Mapping to our Assemblies

When we assemble our genome we lose quantitative information. We can map our reads back to the assembly to determine the relative abundance of each contig in the sample. This is important because we can use this information to bin our sequences into metagenome-assembled genomes (MAGs) - sequences that come from the same organisms should be present in roughly equal proportions. 

We will use two tools - **BWA** and **SAMTools** - to read map our samples. BWA will map short-reads to the assembly. After we will sort the BAM file and index it for rapid random access. 

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

CONTIGS=/home/userid/nematode-microbiome/results/07-assembly
READS=/home/userid/nematode-microbiome/results/03-trimmomatic
OUTPUT=/scratch/userid/nematode-microbiome/results/10-read-map-short-reads

mkdir -p ${OUTPUT}

for FILE in ${CONTIGS}/*; do
  SAMPLE=$(basename ${FILE})
  bwa index ${CONTIGS}/${SAMPLE}/final.contigs.fa
  bwa mem -t 48 ${CONTIGS}/${SAMPLE}/final.contigs.fa ${READS}/${SAMPLE}_R1_paired.fastq.gz ${READS}/${SAMPLE}_R2_paired.fastq.gz | samtools sort -o ${OUTPUT}/alignment.bam --threads 48
  samtools index -@ 48 ${OUTPUT}/alignment.bam
done
```

## Metagenome-assembled Genomes

We are now able to use our abundance information to bin bacterial contigs into metagenome-assembled genomes. We are going to use three tools 1) metabat2, 2) comebin, and 3) dastool. Both Metabat2 and Comebin use tetranucleotide frequencies in conjunction with abundance information for genome reconstruction. However, Comebin uses a contrastive multi-view representation learning to determine the best MAGs and incorporates single-copy gene information and contig length. Finally, Dastool compares the MAGs produced by any binning algorithm and determine the most complete MAGs (with least amount of contamination). 

