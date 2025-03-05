## Taxonomic Composition of Metagenomic Data

There are several tools you can use to assign taxonomy to short metagenomic reads, inlcuding KRAKEN, METPHLAN, and CENETRIFUGE. Today, we will use METAPHLAN4 to assign taxonomy to the nematode microbiome.

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
