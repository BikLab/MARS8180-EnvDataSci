## Metatranscriptomic Data

This week, we will be spending time analyzing a small subset of samples (from 2 stations) that were collected as part of the Tara Oceans Project. There are a total of 12 samples (6 from station TARA_135 and 6 from station TARA_137): 

![TARAOCEANS-CARTE-1024x462](https://github.com/user-attachments/assets/4042a01f-a783-4b37-8c98-97de36bbc751)
![tara-station-map](https://github.com/user-attachments/assets/ad947f27-35fb-47da-919b-23b7f241b157)



1. 3 surface water (5m) samples near Honolulu
2. 3 water samples at the deep-chlorphyll maximum (30-40m) near Honolulu
3. 3 surface water (5m) samples near San Diego
2. 3 water samples at the deep-chlorphyll maximum (30-40m) near San Diego

Example of the metadata file can be seen here: 


| sampleID | sampleENA | event | station | location | long | lat | depth | collection | notes |
|----------|-----------|-------|---------|----------|------|-----|-------|------------|-------|
|ERR1712149| ERS493517	| TARA-20110928Z-SF | TARA_135 | Honolulu | 21.283 | -157.871 | 5m | SEQ-(100L-or-15min)-W>0.8 | surface water | 
| ERR1711927 | ERS493555	| TARA-20110928Z-CH | TARA_135 | Honolulu |21.283 | -157.871 | 30m | SEQ-(100L-or-15min)-W0.8-5	| deep chlorophyll maximum layer |
| ERR1712163	| ERS493652 | TARA-20111124Z-SF | TARA_137 | San Diego | 32.621 | -117.246 | 5m | SEQ-(500mL-or-15min)-N180-2000	| surface water |
| ERR1719158	| ERS493677 | TARA-20111124Z-CH | TARA_137 | San Diego | 32.621 | -117.246 | 40m | SEQ-(100L-or-15min)-W0.8-5 | deep chlorophyll maximum layer |


The steps of analyzing a metatranscriptomic dataset is nearly identical to metagenomic data. We have to 

1. Quality control our raw data
2. Assembly our short-read sequences into contigs
3. Assess quality of our assembly
4. Identify protein coding sequences
5. Assign taxonomy classifcation to each contig
6. Annotate the contigs

**We will largely be following recommendations set forth by Krinos et al (2023) [https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05121-y](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05121-y)**

We have already gone through quality control steps and assembly in the metagenomic section, so we will not spend too much time on that. Instead, we will focus largely on steps 3-6.

## Assess quality of our assembly
We are going to use QUAST to assess the quality of the metatranscriptomic dataset and SALMON to quantiy the expression of transcripts. 

### Running QUAST for Assembly Quality 
Quant is an easy-to-use software that estimates assembly statistics, such as: 

* Number of Contigs
* Length of the Assembly
* N50
* N90
* L50
* L90
* GC%


You can run this tool using the following command: 

```
quast final.contigs.fa -o output/dir --threads 16
```

### Running SALMON for Quantification of Transcripts

Salmon is used to quantify your transcripts. There are two main steps: 

1. index your assembly
2. quantify

These can be run consequetively:

```
salmon index -t final.contigs.fa -i sample-salmon-index -k 31 --threads 12
salmon quant -i sample-salmon-index -l A -1 sample_R1_paired.fastq.gz -2 sample_R2_paired.fastq.gz --validateMappings -o output-sample-dir
```

### Running QUAST and Salmon on Our Data

Let's copy the two scripts over to our home directory

```
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/07-quast.sh /home/userid/metatranscriptomics/scripts
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/08-salmon.sh /home/userid/metatranscriptomics/scripts
```

Let's cd into the script directory to edit and then run these scripts.

```
cd /home/userid/metatranscriptomics/scripts
```

```
nano 07-quast.sh
```

```
#!/bin/sh
#SBATCH --job-name="quast"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 07-quast.err-%N
#SBATCH -o 07-quast.out-%N

module load QUAST

ASSEMBLY=/work/mars8180/instructor_data/metatranscriptome-datasets/results/06-assembly
OUTPUT=/home/userid/metatranscriptomics/results/07-quast

mkdir -p ${OUTPUT}

for folder in ${ASSEMBLY}/*; do
  base=$(basename ${folder})
  quast ${ASSEMBLY}/${base}/final.contigs.fa -o ${OUTPUT}/${base} --threads 16
done
```

The assemblies are located in the instructor_directory. So DO NOT change the file path to the assemblies. Do change the output path to the results folder in your home directory.

Now, we can run Salmon

```
nano 08-salmon.sh
```

```
#!/bin/sh
#SBATCH --job-name="salmon"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 08-salmon.err-%N
#SBATCH -o 08-salmon.out-%N

module load Salmon

ASSEMBLY=/work/mars8180/instructor_data/metatranscriptome-datasets/results/06-assembly
READS=/work/mars8180/instructor_data/metatranscriptome-datasets/results/03-trimmomatic
OUTPUT=/scratch/userid/metatranscriptomics/results/08-salmon

mkdir -p ${OUTPUT}

for folder in ${ASSEMBLY}/*; do
  base=$(basename ${folder})
  salmon index -t ${ASSEMBLY}/${base}/final.contigs.fa -i ${OUTPUT}/${base}/${base}-salmon-index -k 31 --threads 12
  salmon quant -i ${OUTPUT}/${base}/${base}-salmon-index -l A -1 ${READS}/${base}_R1_paired.fastq.gz -2 ${READS}/${base}_R2_paired.fastq.gz --validateMappings -o ${OUTPUT}/${base} --threads 12
done
```

## Identify protein coding sequences

## Assigning taxonomy to each contig

## Annoting the contigs
