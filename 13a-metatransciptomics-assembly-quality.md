## Metatranscriptomic Data

This week, we will be spending time analyzing a small subset of samples (from 2 stations) that were collected as part of the Tara Oceans Project. There are a total of 12 samples (6 from station TARA_135 and 6 from station TARA_137): 

![TARAOCEANS-CARTE-1024x462](https://github.com/user-attachments/assets/4042a01f-a783-4b37-8c98-97de36bbc751)



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

### Assess quality of our assembly
We are going to use QUAST and SALMON to asses the quality of the metatranscriptomic dataset 

### Identify protein coding sequences

### Assigning taxonomy to each contig

### Annoting the contigs
