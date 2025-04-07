## Assesing Quality of Bacterial MAGs

In 2017, two standards were developed by the Genomic Standards Consortium (GSC) for reporting bacterial and archaeal genome sequences. MAGs are classified on different classification levels depending on the completion and contamination determined using single-copy genes.

* Finished Single Contig MAG: >90% complete with less than 5% contamination. Genomes in this category should be a **single contiguous sequence** and also encode the 23S, 16S, and 5S rRNA genes, and tRNAs for at least 18 of the 20 possible amino acids

* High-quality MAG: >90% complete with less than 5% contamination. Genomes in this category should also encode the 23S, 16S, and 5S rRNA genes, and tRNAs for at least 18 of the 20 possible amino acids.

* Medium-quality MAG: ≥50% complete and less than 10% contamination.

* Low-quality MAG: <50% complete with <10% contamination.

It should be noted that there is no minumum assembly size since genomes smaller than 200 kb have been reported. Additionally assembly statistics are usually reported when depositing sequences (N50, L50, largest contig, number of contigs, assembly size, percentage of reads that map back to the assembly, and number of predicted genes per genome).

**Bowers et al. (2017) Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea [https://www.nature.com/articles/nbt.3893](https://www.nature.com/articles/nbt.3893)**

We will use the software program **CheckM2** to determine the quality of the MAGs:

* Chklovski A, Parks DH, Woodcroft BJ, Tyson GW (2023) CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. _Nature Methods_, 20: 1203-1212 - https://www.nature.com/articles/s41592-023-01940-w

It uses two machine learning models to determine the quality depending on wether the genome is **novel** or if it is closely related to a genome in the training set. These two machine learning methods are:  

![41592_2023_1940_Fig1_HTML](https://github.com/user-attachments/assets/ab024524-5f37-40a0-b906-8e2a6f3ad8a2)


1. **Gradient boosted decision trees**: Ke, G. et al. Lightgbm: A highly efficient gradient boosting decision tree. Adv. Neural Inf. Process. Syst. 30, 3146–3154 (2017)
2. **Artifical Neural Networks**: Abadi, M. et al. Tensorflow: a system for large-scale machine learning. In Proc. 12th USENIX Symposium on Operating Systems Design and Implementation (OSDI 16) 265–283 (2016)

The program will determine which machine learning algorithm is best for your genome. 

CheckM2 can be easily run using the following code:

```
checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder>
```



To run CheckM, copy the script from the instructor_directory folder:

```
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts/15-checkm.sh /home/ad14556/nematode-microbiome/scripts/
nano /home/ad14556/nematode-microbiome/scripts/15-checkm.sh
```

```
#!/bin/bash

#SBATCH --job-name="checkm"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e checkm.err-%N
#SBATCH -o checkm.out-%N

# path variables and modules
module load CheckM2

BINS=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/13-dastool-short-reads
OUTPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/15-checkm-short-reads
DATABASE=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/database/CheckM2_database/uniref100.KO.1.dmnd

for FILE in ${BINS}/*; do
  SAMPLE=$(basename ${FILE})
  checkm2 predict --force -x .fa --threads 24 --database_path ${DATABASE} --input ${BINS}/${SAMPLE}/${SAMPLE}_DASTool_bins --output-directory ${OUTPUT}/${SAMPLE}
done
(/home/ad14556
```

## Classifying Bacterial MAGs 

To identify taxonomically classify the the bacterial bins, we can use the GTDB-tk. This tool will identify single-copy genes and place them on a phylogenetic tree. It will out put a text file with the taxonomic id of each bin 

```
cp /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/scripts14-gtdb-tk.sh /home/ad14556/nematode-microbiome/scripts/
nano /home/ad14556/nematode-microbiome/scripts/14-gtdb-tk.sh
```

```
#!/bin/bash

#SBATCH --job-name="gtdb-tk"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e gtdb-tk.err-%N
#SBATCH -o gtdb-tk.out-%N

# path variables and modules
module load GTDB-Tk

BINS=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/13-dastool-short-reads
OUTPUT=/work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/14-gtdbtk-short-reads

for FILE in ${BINS}/*; do
  SAMPLE=$(basename ${FILE})
  gtdbtk classify_wf --genome_dir ${BINS}/${SAMPLE}/${SAMPLE}_DASTool_bins --out_dir ${OUTPUT}/${SAMPLE} --skip_ani_screen -x fa --cpus 24 --pplacer_cpus 24
done
```

