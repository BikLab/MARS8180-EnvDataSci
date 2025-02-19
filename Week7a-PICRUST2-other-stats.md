## PICRUST2 Analysis
PICRUST2 is a method to predict functions using phylogenetic placement of ASV sequences. To do this, we first need to extract the representative sequences and biom table from the QIIME Artifact files. We can do this by using QIIME2 Tools. First, lets copy the script from the `instructor` data directory onto our personal folder. Remember to replace any instances of `userid` with your account ID. 

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
