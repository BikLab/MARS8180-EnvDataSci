## Setting up a projects directory

As discussed earlier in the course, we should first create a directory for our projects. We should have the following subdirectories:

1. data
2. results
3. metadata
4. scripts
5. databases

We can create project directories using the `mkdir` command: 

```
mkdir ddt-project
cd ddt-project
mkdir data databases metadata results scripts
```
Afterwards we are going to copy our data files from the `instructor_data` directory. Login to the transfer node and use the copy command to copy our data files. **DO NOT USE THE MOVE COMMAND**

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/data/* /path/to/your/data/folder
```

## Submitting bash jobs 

Throughout this module, I will show you what commands to use to analyze your dataset. However, you should never run them on a head node. This will slow down the computing node for everyone (and you might receive and email from the GACRC reminding you about good data practices). You should either request an interactive node or submit bash scripts. I recommend the latter for two reasons: 

1. **COMPUTATIONAL RESOURCES**: With bash scripts you can specify the number of memories and CPUs that would help your job run faster
2. **REPRODUCIBILITY**: The scripts in your folder can act as an electronic notebook. Similar to your lab notebook, these functions as notes reminding you what bioinformatics steps you took. 
3. **PUBLISHING**: Is is now common practice to make your code accessible when publishing your manuscript. If you stay organized you can simply copy your scripts folder into a depository (GitHub) and link it in your manuscript's code/data availability subsection.

When you write a bash script you require the following header: 

```
#!/bin/sh 
#SBATCH --job-name="qiime-import"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1:00:00
#SBATCH --mail-user=email@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e file-name.err-%N
#SBATCH -o file-name.out-%N
```

## Demultiplex Data using QIIME2 

## Assessing Data Quality using QIIME2 

First, we need to convert our data into a QIIME2 artifact file (ends with a .qza extension). We can do this using several methods. Either we can use a manifest file that a list of files and their respective paths. Today, we will use the Casava file import. 

**FROM QIIME2**: In Casava 1.8 demultiplexed (paired-end) format, there are two fastq.gz files for each sample in the study, each containing the forward or reverse reads for that sample. The file name includes the sample identifier. The forward and reverse read file names for a single sample might look like `L2S357_15_L001_R1_001.fastq.gz` and `L2S357_15_L001_R2_001.fastq.gz`, respectively. The underscore-separated fields in this file name are:

1. the sample identifier,
2. the barcode sequence or a barcode identifier,
3. the lane number,
4. the direction of the read (i.e. R1 or R2), and
5. the set number.


```
INPUT=/path/to/directory
OUTPUT=/path/to/file/01-qiime-import.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${INPUT} \
  --output-path ${OUTPUT} \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt
```

After we have imported our data as a QZA file type we can use the command `qiime demux summarize` to summarize our data quality. This is going to create a visualization file, indicated by the `.qzv` file extension. 

```
INPUT=/path/to/file/01-qiime-import.qza
OUTPUT=/path/to/file/02-qiime-summarize.qzv

qiime demux summarize \
  --i-data ${INPUT} \
  --o-visualization ${OUTPUT}
```
This command might take a few minutes to run, luckily the instructors thought about this and have the output files saved in the `instructor_data` folder. Copy the output file to your own directory. 

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results/02-qiime-summarize.qzv /path/to/your/file.qza
```
Additionally, we have two more summary files from other sequencing projects that we can use to compare different sequencing quality patterns that you might run into.

**Woodsfall Project:**

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/woodsfall-project/results/02-qiime-summarize.qzv /path/to/your/file.qza
```

**Doliolid Project:**

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/doliolid-project/results/02-qiime-summarize.qzv /path/to/your/file.qza

```

Now we can use the QIIME2 View Tool [https://view.qiime2.org](https://view.qiime2.org) to visualize our data quality. But first, we need to download this to our personal computer. 

**Closely analyze your plots and decide where you would want to trim and truncate your reads and why.**


## Other Tools  
There are multitudes of tools that you can use to visualize and quality control your data. Two of them are fastQC and multiQC. Kevin let us borrow his data (illumina sequnecing of vibrio isolates) to analze in class. 

To run fastqc you can use the following command.

```
module load FastQC

INPUT=/path/to/data/directory
OUTPUT=/path/to/results/directory/01-fastqc

mkdir -p ${OUTPUT}
fastqc ${INPUT}/* -o ${OUTPUT} -t 8
```
The flag `-t` specifies the number of threads. However, make sure their is a minimum of 150Mb allocated for each thread. 

Then we can run multiqc which aggregates all of our fastqc files into one readable HTML file. 

```
module load MultiQC

INPUT=/path/to/results/directory/01-fastqc
OUTPUT=/path/to/results/directory/02-multiqc

multiqc --outdir ${OUTPUT} ${INPUT}
```

Copy the file to your personal computer and open the `muiltiqc_report.html` file
 
