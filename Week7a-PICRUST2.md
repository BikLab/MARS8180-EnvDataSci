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


To get started, we first need to extract the representative sequences and biom table from the QIIME Artifact files. We can do this by using QIIME2 Tools. First, lets copy the script from the `instructor` data directory onto our personal folder. Remember to replace any instances of `userid` with your account ID. 

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/scripts/08-extract-data-qiime2.sh /home/userid/ddt-project/scripts
```

Now we can nano into the folder and edit our paths 

```
nano /home/userid/ddt-project/scripts
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

REFSEQ=/home/userid/ddt-project/results/05-dada2-rep-seq.qza
REFTABLE=/home/userid/ddt-project/results/05-dada2-feature-table.qza
SEQOUT=/home/userid/ddt-project/results/08-repseq/
TABLEOUT=/home/userid/ddt-project/results/08-biomtable/
```
module load QIIME2

qiime tools export \
  --input-path ${REFSEQ} \
  --output-path ${SEQOUT}

qiime tools export \
  --input-path ${REFTABLE} \
  --output-path ${TABLEOUT}
```
make sure that your paths are correct and pointing to the correct locations. The variable `REFSEQ` and the `REFTABLE` should be path to the ASV representative sequences and the feature table output by the DADA2 script. 

Now we can submit the job to the cluster. This will generate two folders with our exported sequences (`08-repseq/`) and the biom table (`08-biomtable`). We can use these files to run the picrust2 script and predict functions based on our phylogenetic placement of the ASVs.
