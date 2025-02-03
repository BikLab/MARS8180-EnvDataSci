### Week 5a Clustering eDNA metabarcoding data

In today's class we'll cover the next two steps in the bioinformatics pipeline
1. Adaptor Trimming & Read Merging
2. Clustering raw eDNA reads into Amplicon Sequence Variants (ASVs)

#### 11:10 Minute Cards

Fill out Minute cards: https://forms.gle/fK2FGG1uUSoaZTSo6

#### 11:10-11:25: Adaptor Trimming

First let's log into the teaching cluster and then copy the scripts for adapter trimming and denoising using DADA2

```
ssh ugaid@teach.gacrc.uga.edu
```

```
cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/scripts/03-cutadapt.sh /home/ad14556/ddt-project/scripts/

cp /work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/scripts/03-denoise-dada2.sh /home/ad14556/ddt-project/scripts/

```

We will use a tool called cutadapt to remove our primers and adapters. After we will summarize the sequences that were successfully trimmed. We will request 12 CPUs to multithread cutadapt. 

```
INPUT=/home/ugaid/ddt-project/results/01-qiime-import.qza
OUTPUT=/home/ugaid/ddt-project/results/03-cutadapt-sans-primers.qza
SUMM=/home/ugaid/ddt-project/results/03-cutadapt-summ.qzv

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ${INPUT} \
  --p-adapter-f GTGYCAGCMGCCGCGGTAA \
  --p-adapter-r GGACTACNVGGGTWTCTAAT \
  --p-error-rate 0.1 \
  --o-trimmed-sequences ${OUTPUT} \
  --p-cores 12 \
  --verbose
  
qiime demux summarize \
  --i-data ${OUTPUT} \
  --o-visualization ${SUMM}
```

#### 11:25-11:45 Clustering eDNA reads into Amplicon Sequence Variants

Graphics taken from Ben Callahan's lectures at the STAMPS 2024 course at MBL: https://github.com/mblstamps/stamps2024/wiki#17

Different names for eDNA "clusters", which are all essentially a molecular proxy for a biological species (see also slide further down)
* OTU = Operational Taxonomic Unit
* ASV = Amplicon Sequence Variant
* ESV = Exact Sequence Variant
* zOTU = Zero-radius OTU (a denoised sequence)

---
<img width="1151" alt="Screenshot 2025-02-02 at 2 12 53 PM" src="https://github.com/user-attachments/assets/cbee35f6-cb9a-42ce-a0ea-82e8b2a5da82" />

---
<img width="1153" alt="Screenshot 2025-02-02 at 2 12 19 PM" src="https://github.com/user-attachments/assets/74fcd062-7725-4d66-8d8c-bd8e702361c9" />

---
<img width="1125" alt="Screenshot 2025-02-02 at 2 13 45 PM" src="https://github.com/user-attachments/assets/33dcb315-aada-48a0-ade9-a36e9cfffb88" />

---
Previously, we used to cluster eDNA reads under some % identity cutoff (e.g. 97% Operational Taxonomic Units), by either
1) comparing against a database and discarding reads that didn't match the database (closed-reference OTUs)
2) sequential database comparison followed by de novo clustering of non-database matches (open-reference OTUs)
3) De novo clustering of all sequences by comparing eDNA reads against each other (de novo OTUs)

The nucleotide sequence chosen to represent all reads in an OTU is known as a "representative sequence"

<img width="1137" alt="Screenshot 2025-02-02 at 2 15 09 PM" src="https://github.com/user-attachments/assets/2f6e06ab-32a9-4ed8-85a3-db5b1cd529c8" />

---
<img width="1149" alt="Screenshot 2025-02-02 at 2 21 57 PM" src="https://github.com/user-attachments/assets/4b243148-6d9d-49dd-8c83-d985cc0c8c1f" />

McLaren & Callahan, mBio 2018 - https://journals.asm.org/doi/10.1128/mbio.02149-17 

---
<img width="1149" alt="Screenshot 2025-02-02 at 2 23 10 PM" src="https://github.com/user-attachments/assets/d1c142d0-7fd5-4325-afec-61569894e185" />

Callahan, McMurdie & Holmes 2017: https://www.nature.com/articles/ismej2017119

---
<img width="1144" alt="Screenshot 2025-02-02 at 2 24 15 PM" src="https://github.com/user-attachments/assets/d6875e8f-28cc-472d-851d-572d02f8b28e" />

---
<img width="1151" alt="Screenshot 2025-02-02 at 2 25 04 PM" src="https://github.com/user-attachments/assets/441960ab-63df-457d-9c97-bd38de830b3d" />

---
<img width="1147" alt="Screenshot 2025-02-02 at 2 27 37 PM" src="https://github.com/user-attachments/assets/d7bee1a6-8144-4651-85f8-357cc24a08d3" />

---
<img width="1152" alt="Screenshot 2025-02-02 at 2 27 56 PM" src="https://github.com/user-attachments/assets/851ae1a1-3697-4b83-a764-86c71bca1e97" />

---
<img width="1150" alt="Screenshot 2025-02-02 at 2 28 14 PM" src="https://github.com/user-attachments/assets/5488787a-9856-4839-a436-e64d56f32ae2" />

---
<img width="1144" alt="Screenshot 2025-02-02 at 2 29 17 PM" src="https://github.com/user-attachments/assets/e37459ad-734d-4c65-b474-1d2e68f359ab" />

---
Ben Callahan's reccomendations: 
<img width="1152" alt="Screenshot 2025-02-02 at 2 30 05 PM" src="https://github.com/user-attachments/assets/cc6d849f-9d98-4fe9-9bcb-d68cbbf3fbed" />

---
#### 11:45-11:50 Class Reflection

What issues or biological uncertainties arise when using eDNA reads and ASVs instead of "traditional approaches" (microscopy, cell counts, net collections, etc.)? Take 5 minutes to jot down some thoughts here: https://docs.google.com/document/d/1vUWVESKF8idcMElPR1QVl3EByDOU-6UmsD55VNdKN_s/edit?usp=sharing

#### 11:50 - 12:25 Clustering eDNA data using DADA2

We will denoise our dataset using the DADA2 software using 12 threads.

```
INPUT=/home/ad14556/ddt-project/results/03-cutadapt-sans-primers.qza
OUTPUT=/home/ad14556/ddt-project/results

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${INPUT} \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences ${OUTPUT}/03-dada2-rep-seq.qza \
  --o-table ${OUTPUT}/03-dada2-feature-table.qza \
  --o-denoising-stats ${OUTPUT}/03-dada2-stats.qza \
  --p-n-threads 12
```
