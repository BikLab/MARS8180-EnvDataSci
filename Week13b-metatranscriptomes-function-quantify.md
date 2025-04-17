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
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/11-transdecoder.sh /home/userid/metatranscriptomics/scripts
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

ASSEMBLY=/work/mars8180/instructor_data/metatranscriptome-datasets/results/09-merged-assembly/merged-assembly-98-simplify-namese.fa
OUTPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/11-transdecoder

mkdir -p ${OUTPUT}

TransDecoder.LongOrfs -t ${ASSEMBLY} -m 100 --output_dir ${OUTPUT}
TransDecoder.Predict -t ${ASSEMBLY} --output_dir ${OUTPUT} --no_refine_starts
```

## Rerun salmon to quantify our data

Before we rerun Salmon, we need to fix the contig/sequence names in our files so that that are 1) simplified and 2) are unique. We can do this by using a couple of bash tools: 

```
#remove white spaces in contig names
sed 's/ /_/g' merged-assembly-98.fa > merged-assembly-98-sans-space.fa 

# simplify the contig names so that they start with contig and are. sequential
seqtk rename merged-assembly-98-sans-space.fa contig_ > merged-assembly-98-simplify-namese.fa
```

Then we can run Salmon; however, there are a few differences. This time we only need to index one file (the merged assembly). AND we need to merge the quant.sf files together so that we can directly compare samples.

Lets copy the script to our home directory

```
cp /work/mars8180/instructor_data/metatranscriptome-datasets/scripts/10-salmon.sh /home/userid/metatranscriptomics/scripts
```

```
nano /home/userid/metatranscriptomics/scripts/10-salmon.sh
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
#SBATCH -e 10-salmon.err-%N
#SBATCH -o 10-salmon.out-%N

module load Salmon
module load seqtk

ASSEMBLY=/work/mars8180/instructor_data/metatranscriptome-datasets/results/09-merged-assembly
READS=/work/mars8180/instructor_data/metatranscriptome-datasets/results/03-trimmomatic
OUTPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/10-salmon

mkdir -p ${OUTPUT}

sed 's/ /_/g' ${ASSEMBLY}/merged-assembly-98.fa > ${ASSEMBLY}/merged-assembly-98-sans-space.fa #remove white spaces in contig names
seqtk rename ${ASSEMBLY}/merged-assembly-98-sans-space.fa contig_ > ${ASSEMBLY}/merged-assembly-98-simplify-namese.fa # simplify the contig names

salmon index -t ${ASSEMBLY}/merged-assembly-98-simplify-namese.fa -i ${OUTPUT}/merged-salmon-index -k 31 --threads 12

for file in ${READS}/*_R1_paired.fastq.gz; do
  base=$(basename ${file} _R1_paired.fastq.gz)
  salmon quant -i ${OUTPUT}/merged-salmon-index -l A -1 ${READS}/${base}_R1_paired.fastq.gz -2 ${READS}/${base}_R2_paired.fastq.gz --validateMappings -o ${OUTPUT}/${base} --threads 12
done

salmon quantmerge --quants ${OUTPUT}/ERR* -o ${OUTPUT}/salmon-all-sample-quant.txt
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

<img width="645" alt="Screenshot 2025-04-10 at 10 20 24 AM" src="https://github.com/user-attachments/assets/882eb0b8-6d69-4a36-92a9-4380ea105fe2" />




First let's install this using a combination of conda and pip

```
interact --mem=15G

mkdir /home/ad14556/conda-env/eukulele
conda create -p /home/ad14556/conda-env/eukulele
source activate /home/ad14556/conda-env/eukulele

pip install EUKulele
```

Now we are able to run Eukulele on our metatranscriptomes

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

module load EUKulele
module load DIAMOND

INPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/scripts
OUTPUT=/work/mars8180/instructor_data/metatranscriptome-datasets/results/12-eukulele
DATABASE=/work/mars8180/instructor_data/metatranscriptome-datasets/databases

mkdir -p ${OUTPUT}

EUKulele --mets_or_mags mets -s ${INPUT} --p_ext ".pep" -d phylodb --reference_dir ${DATABASE} -o ${OUTPUT} --CPUs 12
```


## Annoting the contigs
