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

Now we can download the html files on our personal computers to assess transcript assembly quality. 

## Merging Transcriptomes


```
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/09-merge-assemblies.sh /home/userid/metatranscriptomics/scripts/
```

```
nano /home/userid/metatranscriptomics/scripts/09-merge-assemblies.sh
```

```
#!/bin/sh
#SBATCH --job-name="merge"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 09-merge.err-%N
#SBATCH -o 09-merge.out-%N

module load MMseqs2

INPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/06-assembly
OUTPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/09-merged-assembly

mkdir -p ${OUTPUT}
cat ${INPUT}/*/*.fa > ${OUTPUT}/merged-assembly.fa

mmseqs createdb ${OUTPUT}/merged-assembly.fa ${OUTPUT}/merged-assembly
mmseqs linclust ${OUTPUT}/merged-assembly ${OUTPUT}/merged-assembly-2 tmp --min-seq-id 0.98 --cov-mode 1 --split-memory-limit 120G --remove-tmp-files
mmseqs createsubdb ${OUTPUT}/merged-assembly-2 ${OUTPUT}/merged-assembly ${OUTPUT}/merged-assembly-3
mmseqs convert2fasta ${OUTPUT}/merged-assembly-3 ${OUTPUT}/merged-assembley-98.fa
```



## Identify protein coding sequences
To identify protein coding regions within transcripts we will use the tool **TransDecoder**. TransDecoder identifies likely coding sequences based on the following criteria:

1. a minimum length open reading frame (ORF) is found in a transcript sequence
2. a log-likelihood score similar to what is computed by the GeneID software is > 0.
3. the above coding score is greatest when the ORF is scored in the 1st reading frame as compared to scores in the other 2 forward reading frames.
4. if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF, the longer one is reported. However, a single transcript can report multiple ORFs (allowing for operons, chimeras, etc).

TransDecoder is run in two steps. First, sequences with long open reading frames are extracted from the assembled sequences.  By default, only sequences that are 100 amino acids long are kept. Second, you predict the likely coding region. 

In practice, it looks like the following:

```
TransDecoder.LongOrfs -t final.contigs.fa -m 100 --output_dir output/dir/sample-basename
TransDecoder.Predict -t final.contigs.fa --output_dir output/dir/sample-basename
```

You will have four final outputs: 

* **transcripts.fasta.transdecoder.pep**: peptide sequences for the final candidate ORFs; all shorter candidates within longer ORFs were removed. This is the most important file you'll need.

* **transcripts.fasta.transdecoder.cds**: nucleotide sequences for coding regions of the final candidate ORFs
 
* **transcripts.fasta.transdecoder.gff3**: positions within the target transcripts of the final selected ORFs

* **transcripts.fasta.transdecoder.bed**: bed-formatted file describing ORF positions, best for viewing using GenomeView or IGV.


### Running TransDecoder on our data

First, lets copy the script to our home directory

```
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/09-transdecoder.sh /home/userid/metatranscriptomics/scripts
```

```
nano /home/userid/metatranscriptomics/scripts/11-transdecoder.sh
```

```
#!/bin/sh
#SBATCH --job-name="transdecoder"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 09-transdecoder.err-%N
#SBATCH -o 09-transdecoder.out-%N

module load TransDecoder

ASSEMBLY=/work/mars8180/instructor_data/metatranscriptome-datasets/results/06-assembly
OUTPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/09-transdecoder

mkdir -p ${OUTPUT}

for folder in ${ASSEMBLY}/*; do
  base=$(basename ${folder})
  mkdir -p ${OUTPUT}/${base}
  cd ${OUTPUT}/${base}
  TransDecoder.LongOrfs -t ${ASSEMBLY}/${base}/final.contigs.fa -m 100 --output_dir ${OUTPUT}/${base}
  TransDecoder.Predict -t ${ASSEMBLY}/${base}/final.contigs.fa --output_dir ${OUTPUT}/${base} --no_refine_starts
done
```


## Assigning taxonomy to each contig
We will be using Eukulele and the PhyloDB database [https://github.com/allenlab/PhyloDB](https://github.com/allenlab/PhyloDB) to assign taxonomy to our contig. 

First, lets install this using conda and pip

```
interact --mem=15G

mkdir /home/userid/conda-env/eukulele
conda create -p /home/userid/conda-env/eukulele
source activate /home/userid/conda-env/eukulele
pip install EUKulele
```

I have already installed the phyloDB database in the instructory directory `/work/mars8180/instructor_data/metatranscriptome-datasets/databases/phylodb` using the following command

```
EUKulele download --database phylodb
```

Now, we can classify taxonomy using EUKulele

```
EUKulele --mets_or_mags mets --sample_dir ${TRANSDECODER-DIRECTORY} --p_ext ".pep" --database ${DATABASE} -o sample/eukulele-output
```

### Running EUKulele on our data

Let's run this on our metatranscriptomes

```
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/12-eukulele.sh /home/userid/metatranscriptomics/scripts

nano /home/userid/metatranscriptomics/scripts/12-eukulele.sh
```

```
#!/bin/sh
#SBATCH --job-name="eukulele"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 10-eukulele.err-%N
#SBATCH -o 10-eukulele.out-%N

module load Miniconda3
source activate /home/ad14556/conda-env/eukulele

INPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/09-transdecoder
OUTPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/10-eukulele
DATABASE=/work/mars8180/instructor_data/metatranscriptome-datasets/databases/phylodb

mkdir -p ${OUTPUT}

for folder in ${INPUT}/*; do
  base=$(basename ${folder})
  EUKulele --mets_or_mags mets --sample_dir ${INPUT}/${base} --p_ext ".pep" --database ${DATABASE} -o ${OUTPUT}/${base}
done
```


## Annoting the contigs
