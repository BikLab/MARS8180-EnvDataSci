## Setting up a CONDA environment
First, lets request an interactive node with 12 GB of memory as recommended by CONDA.

```
srun --pty  --cpus-per-task=1 --job-name=interact \
--ntasks=1 --nodes=1 --partition=batch --time=02:00:00 --mem=12GB /bin/bash -l
```

Now we can load the Miniconda module and create a conda repository in your home directory.

```
module load Miniconda3
mkdir conda-env
```

Then, lets create a repository for our qiime2 conda environment. The name should contain the version of the environment we are trying to build. 

```
mkdir conda-env/qiime2-amplicon-2024.10
conda create -p conda-env/qiime2-amplicon-2024.10
```

We should download the conda yaml file using the curl command. We are going to use a new parameter we have not seen before `-L`. This allows any redirects from the URL. 

```
curl -LO "https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml"
```

Now we can load the environment and install qiime2 

```
source activate /home/ad14556/conda-env/qiime2-amplicon-2024.10
conda env update --file qiime2-amplicon-2024.10-py310-linux-conda.yml
```

While your environment is being set up, let's navigate to the anaconda website [https://anaconda.org](https://anaconda.org). We can search for packages that are distributed using the conda package management system. Search for the QIIME2 software. 

How many variations of QIIME2 are there? Why would we want to create a conda environment for specifically for the QIIME2 amplicon package? 

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
Afterwards we are going to copy our data files from the `instructor_data` directory. Login to the transfer node and use the copy command to copy our data files and use the `tar` command to unzip and decompress our file. **DO NOT USE THE MOVE COMMAND**

```
cd data
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/data/ddt-raw-fastq.tar.gz .
tar -xvzf ddt-raw-fastq.tar.gz
```

## Submitting bash jobs 

Throughout this module, I will show you what commands to use to analyze your dataset. However, you should never run them on a head node. This will slow down the computing node for everyone (and you might receive and email from the GACRC reminding you about good data practices). You should either request an interactive node or submit bash scripts. I recommend the latter for two reasons: 

1. **COMPUTATIONAL RESOURCES**: With bash scripts you can specify the number of memories and CPUs that would help your job run faster
2. **REPRODUCIBILITY**: The scripts in your folder can act as an electronic notebook. Similar to your lab notebook, these scripts are notes reminding you what bioinformatics steps you took. 
3. **PUBLISHING**: It is now common practice to make your code accessible when publishing your manuscript. If you stay organized you can simply copy your scripts folder into a depository (GitHub) and link it in your manuscript's code/data availability subsection.

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
The sequencing facility will give you your data as 1) multiplexed or 2) demultiplexed sequences. If your sequencing data is a multiplexed file, all of your sample sequencing data will be located in a single file. Demultiplexing is a process of seperating your samples using your unique sequencing barcodes. 

**If the sequencing facility demultiplexes your sequencing data, they will also send you the multiplexed raw sequencing data**

If you have multiplexed paired-end samples with barcodes in separate fastq file, you will typically recieve 4 files: 

1. Forward sequences (`project_L001_R1.fastq.gz`)
2. Reverse sequences (`project_L001_R2.fastq.gz`)
3. Forward barcodes (`barcodes_L001_I1.fastq.gz`)
4. Reverse barcodes (`barcodes_L001_I2.fastq.gz`)

Additionally, you can have multiplexed paired-end samples with barcodes within your sequences. In this case, you will have two sequencing files: 

1. Forward sequences (`project_L001_R1.fastq.gz`)
2. Reverse sequences (`project_L001_R2.fastq.gz`)

If you are working with single-end data you will receive the following files: 

1. A single sequences (`project_L001.fastq.gz`)
2. A single barcodes file (`barcodes.fastq.gz`) 

Regardless of the format, you will need to use the barcodes in the metadata file to indicate which barcodes are associated with which samples. **Note: Nowadays, the sequencing facility typically demultiplexes your sample free of cost. However, you must share your metadata and barcodes with the facility**

A typically metadata file with barcode(s) for paired-end sequence looks like this: 

 
|SampleID|BarcodeSequence|LinkerPrimerSequence|ReverseBarcode|
|--------|---------------|--------------------|--------------|
|sampleA |GCTAGCCTTCGTCGC|TATGGTAATTGTGTGYCAGCMGCCGCGGTAA|GATCGGGACACCCGA|
|sampleB |GCTAGCCTTCGTCGC|TATGGTAATTGTGTGYCAGCMGCCGCGGTAA|GATCTGTCTATACTA|
|sampleC |GCTCCTAACGGTCCA|TATGGTAATTGTGTGYCAGCMGCCGCGGTAA|GATAATAACTAGGGT|
|sampleD |GCTCGCGCCTTAAAC|TATGGTAATTGTGTGYCAGCMGCCGCGGTAA|GATTACGGATTATGG|
|sampleE |GCTACTACTGAGGAT|TATGGTAATTGTGTGYCAGCMGCCGCGGTAA|GATGTGGAGTCTCAT|


## What is a FASTQ file
A fastq file is a text file that stores information about the sequence and its quality score of each base. 

A typical fastq file looks like the following: 

```
@SRR2584863.1 HWI-ST957:244:H73TDADXX:1:1101:4712:2181/1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGCAGGGGTTTTTGGGTTACCTGATCCTGAGAGTTAACGGTAGAAACGGTCAGTACGTCAGAATTTACGCGTTGTTCGAACATAGTTCTG
+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################
```

1. The first line always starts with `@` and then information about the read - the instrument name, the coordinates of the flow cells, and the members of a pair. Typically, `/1` refers to a forward read and `/2` are the reverse reads. 
2. The second line is the raw sequences.
3. The third line is `+` character that seperates the sequences and sequence quality score.
4. The last line, is the sequence quality encoded by characters that represent a specific PHRED score.

From lowest (99.999999% probability of an error) to highest quality score (0.0001% probability of an error): 
```!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI```

Lets use the following command to get the first 8 lines of one of our ddt-samples. 

```
zcat ddt-raw-fastq/DDT.10.1_S1258_L001_R1_001.fastq.gz | head -n 8
```

```
@VH00301:398:AAFT7FLM5:1:1101:23893:1095 1:N:0:CTTCAAGATTTC+ATGATACGTAAT
ATGCACGTCCCAGTATTGGTTTCAACATTGTGAAACTGTGAATTGCTCAATAAAACAGCTATTGCATATATGACTGCATATTTACATGGCTAGCCGTGGCAATTCTAGAGCTAATACATGTACGGAGCCTAACTTTGTGGGGAGGGTATTGTTTATTAGTTGTGGAACCAGTCCAGGTTATCCTTGGTTTCTTCTGATTGGTTGTAACTAAATGAATCGCATGGCATCAGTTGGTGGTGCATCATTCAAGCATCTGACCTATCAGCTTCCGACGGTAAGGTATTGGCCTACCTTGGCAATG
+
C5CCCCCCCCCCC*CCCCCCCCCCCCCCCCCCCCCCCC*CCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCC*CCCCCCCCCCCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCCC5CCCCCCCCCCCCCC
@VH00301:398:AAFT7FLM5:1:1101:39761:1170 1:N:0:CTTCAAGATTTC+ATGATACGTAAT
GTGCATGTCTCAGTATAAGTGTTTCACTGCGAAACTGCGAATGGCTCATTAAAACAGTTATAGTTTCCATGTCAGTTGTTTATTACCTGGATATCCACGGTAATTCTAGAGCTAATACATGCGTCCAAACCCGACTTTTGCGGAAGGGTTGTGCTTATTAGACACTGAACCATCCCGGGCTTGCCCGGTTTCGAGGTGATTGATGGTAATCGAACGAATCGCATGCTTCGGCGGCGATGATTCATTCAAGTTTCTGACCTATCAGCTTCCGACGGTAGGGTATTGGCCTACCGTGGCTTTG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC5CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*
```


## Assessing Data Quality using QIIME2 

First, we need to convert our data into a QIIME2 artifact file (ends with a .qza extension). We can do this using several methods. For example, we can use a manifest file that a list of files and their respective paths. 

A manifest will contain a sample identifier, along with the file paths of those samples. The manifest file can contain variables (e.g., $HOME or $PWD). The following example illustrates a simple fastq manifest file for paired-end read data for three samples.


|sample-id|forward-absolute-filepath|reverse-absolute-filepath|
|---------|-------------------------|-------------------------|
|sampleA|$PWD/some/filepath/sampleA_R1.fastq.gz|$PWD/some/filepath/sampleA_R2.fastq.gz|
|sampleB|$PWD/some/filepath/sampleB_R1.fastq.gz|$PWD/some/filepath/sampleB_R2.fastq.gz|
|sampleC|$PWD/some/filepath/sampleB_R1.fastq.gz|$PWD/some/filepath/sampleC_R2.fastq.gz|


However, today, we will use the Casava file import. 

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
This command might take a few minutes to run, luckily the instructors thought about this and have the output files saved in the `instructor_data` folder. Copy the output file to your own computer. 

```
scp userid@txfer.gacrc.uga.edu:"/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/results/02-qiime-summarize.qzv" /file/path/02-qiime-summarize-ddt.qzv
```
Additionally, we have two more summary files from other sequencing projects that we can use to compare different sequencing quality patterns that you might run into.

**Woodfall Project:**

```
scp userid@txfer.gacrc.uga.edu:"/work/mars8180/instructor_data/metabarcoding-datasets/woodsfall-project/results/02-qiime-summarize.qzv" /file/path/02-qiime-summarize-woodfall.qzv
```

**Doliolid Project:**

```
scp userid@txfer.gacrc.uga.edu:"/work/mars8180/instructor_data/metabarcoding-datasets/doliolid-project/results/02-qiime-summarize.qzv" /file/path/02-qiime-summarize-doliolid.qzv

```

Now we can use the QIIME2 View Tool [https://view.qiime2.org](https://view.qiime2.org) to inspect and visualize our sequencing data quality. But first, we need to download our `.qzv` data files to our personal computer. 

**Closely analyze your plots and decide where you would want to trim and truncate your reads and why.**


## Other Tools  
There are multitudes of tools that you can use to visualize and quality control your data. Two of them are fastQC and multiQC. Kevin let us borrow his data (Illumina sequencing of vibrio isolates) to analyze them in class. 

To run fastqc you can use the following command.

```
module load FastQC

INPUT=/path/to/data/directory
OUTPUT=/path/to/results/directory/01-fastqc

mkdir -p ${OUTPUT}
fastqc ${INPUT}/* -o ${OUTPUT} -t 8
```
The flag `-t` specifies the number of threads. However, we need to make sure we allocate a minimum of 250Mb to each thread. 

Then we can run multiqc which aggregates all of our fastqc files into one readable HTML file. 

```
module load MultiQC

INPUT=/path/to/results/directory/01-fastqc
OUTPUT=/path/to/results/directory/02-multiqc

multiqc --outdir ${OUTPUT} ${INPUT}
```

Copy the file to your personal computer and open the `muiltiqc_report.html` file
 
