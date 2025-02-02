### Week 5a Clustering eDNA metabarcoding data

In today's class we'll cover the next two steps in the bioinformatics pipeline
1. Adaptor Trimming & Read Merging
2. Clustering raw eDNA reads into Amplicon Sequence Variants (ASVs)

#### 11:10 Minute Cards

Fill out Minute cards: https://forms.gle/fK2FGG1uUSoaZTSo6

#### 11:10-11:25: Adaptor Trimming & Read Merging

ALEJANDRO INSERT CONTENT + CODE HERE

#### 11:25-11:45 Clustering eDNA reads into Amplicon Sequence Variants

Graphics taken from Ben Callahan's lectures at the STAMPS 2024 course at MBL: https://github.com/mblstamps/stamps2024/wiki#17

Different names for eDNA "clusters", which are all essentially a molecular proxy for a biological species:
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

<img width="1137" alt="Screenshot 2025-02-02 at 2 15 09 PM" src="https://github.com/user-attachments/assets/2f6e06ab-32a9-4ed8-85a3-db5b1cd529c8" />


#### 11:45-11:50 Class Reflection

What issues or biological uncertainties arise when using eDNA reads and ASVs instead of "traditional approaches" (microscopy, cell counts, net collections, etc.)? Take 5 minutes to jot down some thoughts here: https://docs.google.com/document/d/1vUWVESKF8idcMElPR1QVl3EByDOU-6UmsD55VNdKN_s/edit?usp=sharing

#### 11:50 - 12:25 Clustering eDNA data using DADA2

ALEJANDRO INSERT CONTENT + CODE HERE

