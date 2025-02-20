## PICRUST2 Analysis
PICRUST2 is a method to predict functions using phylogenetic placement of ASV sequences - it can help you explore the **functional potential** of a microbial community, which can then (ideally) be tested/validated through other types of -Omics sequencing such as metagenomics or metatranscriptomics. Note: Some people are skeptical of PICRUSt2 results - as they are highly dependent on reference genome databases - but we still think it is a useful (and popular) tool. 
* Github wiki with tutorials - https://github.com/picrust/picrust2/wiki
* Paper to cite: Douglas, G.M., Maffei, V.J., Zaneveld, J.R. et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38, 685–688 (2020). https://doi.org/10.1038/s41587-020-0548-6

![PICRUSt2_flowchart](https://github.com/user-attachments/assets/e072d959-7c83-48e1-848c-607646cc45da)
(Figure from PICRUSt2 wiki)

---

There are several steps to the PICRUST workflow: 

1. Align ASVs to Reference Sequences
2. Place ASVs to Reference Trees
3. Infer Copy Number
4. Metagenome Predications Based on Phylogenetic Placement
5. Infer Pathway Abundances

However, we can easily run PICRUST2 using a single command that wraps all of our steps into one workflow. 

You can only run PICRUST on 16S rRNA metabarcoding datasets, so to show you how to run this pipeline we are going to copy over a projects folder from the instructor data folder. Lets copy the script from the `instructor` data directory onto our personal folder. Remember to replace any instances of `userid` with your account ID. 

```
cp -r /work/mars8180/instructor_data/metabarcoding-datasets/memb-project home/userid/
```

Now we can nano into the scripts and edit our paths. First, lets nano into the script names `07-extract-bio-table-ref-seq.sh`. This script will allow us to export the data (representative sequences and biom table from our QIIME2 artifact files. 

```
nano /home/userid/ddt-project/scripts/07-extract-bio-table-ref-seq.sh
```

```
#!/bin/sh
#SBATCH --job-name="extract-data"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 07-extract-data.err-%N
#SBATCH -o 07-extract-data.out-%N

module load QIIME2/2024.10-amplicon

REFSEQ=/home/userid/memb/results/dada2-rep-seqs.qza
REFTABLE=/home/userid/memb/results/dada2-table.qza
SEQOUT=/home/userid/memb/results/repseq/
TABLEOUT=/home/userid/memb/results/biomtable/

module load QIIME2

qiime tools export \
  --input-path ${REFSEQ} \
  --output-path ${SEQOUT}

qiime tools export \
  --input-path ${REFTABLE} \
  --output-path ${TABLEOUT}
```

Make sure that your paths are pointing to the correct locations. The variable `REFSEQ` and the `REFTABLE` should be path to the ASV representative sequences and the feature table output by the DADA2 script. 

Now we can submit the job to the cluster. This will generate two folders with our exported sequences (`08-repseq/`) and the biom table (`08-biomtable`). We can use these files to run the picrust2 script and predict functions based on our phylogenetic placement of the ASVs.

This time for this scripts, we are going to use 12 threads and 10Gb per thread.

```
nano /home/userid/ddt-project/scripts/08-picrust2.sh
```

```
#!/bin/bash
#SBATCH --job-name="08-picrust2"
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e 08-picrust2.err-%N
#SBATCH -o 08-picrust2.out-%N

#Path Variables
SEQ=/home/ad14556/memb/results/repseq/dna-sequences-rc.fasta
BIOM=/home/ad14556/memb/results/biomtable/feature-table.biom
OUT=/home/ad14556/memb/results/08-picrust/

# Activate Picrust2
module load PICRUSt2

# Picrust2 commands
picrust2_pipeline.py \
  -s ${SEQ} \
  -i ${BIOM} \
  -o ${OUT} -p 12 \
  --stratified \
  --in_traits EC,KO
```

The key output files are:

1. EC_metagenome_out - Folder containing unstratified EC number metagenome predictions (pred_metagenome_unstrat.tsv.gz), sequence table normalized by predicted 16S copy number abundances (seqtab_norm.tsv.gz), and the per-sample NSTI values weighted by the abundance of each ASV (weighted_nsti.tsv.gz).
2. KO_metagenome_out - As EC_metagenome_out above, but for KO metagenomes.
3. pathways_out - Folder containing predicted pathway abundances and coverages per-sample, based on predicted EC number abundances.
