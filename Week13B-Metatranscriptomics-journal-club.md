## Metatranscriptomics Journal Club

_"A **metatranscriptome** represents the collection of all RNA transcripts (e.g., mRNA, rRNA, tRNA) produced by all the microorganisms in a particular environment at a given point in time, and provides insights into the gene expression patterns and activities of the microbial community."_ (Aplakidou et al. 2024)

Moran MA, Satinsky B, Gifford SM, Luo H, Rivers A, Chan LK, Meng J, Durham BP, Shen C, Varaljay VA, Smith CB (2013) Sizing up metatranscriptomics. _ISME_, 7(2):237-43. https://www.nature.com/articles/ismej201294

Metatranscriptomics studies are inherently temporal in nature: _"A typical marine bacterial cell in coastal seawater contains only ~200 molecules of mRNA, each of which lasts only a few minutes before being degraded. Such a surprisingly small and dynamic cellular mRNA resevoir has important implications for understanding the bacterium's responses to environmental signals, as well as for our ability to measure those responses."_ (Abstract of Moran et al. 2013)

A visual on what this looks like in the environment from this paper: 

<img width="326" alt="Screenshot 2025-04-02 at 6 16 31 PM" src="https://github.com/user-attachments/assets/778b72ba-5839-4e2b-b705-8d2e41a3d5f4" />

* Lab studies indicate that **mRNA has an average half-life of 2.4 to 5 minutes**. In contrast, the half-life of a typical bacterial protein is ~20 hours, which is about two orders of magnitude longer than an mRNA half-life.
* mRNA responds sensitively to the beginning and end of envionrmental fluctions - relative shifts in proteins show a lag time, and proteins are thus less responsive to rapid environmental dynamics.

#### Four Key Observations from Moran et al. 2013: 

1. The abundance of mRNAs from functional genes is not a reliable rate proxy for those functions in naturally fluctuating environments, and neither is the abundance of proteins
2. Instantaneous inventories of mRNA pools are nontheless highly informative about ongoing ecologically relevant processes
3. Fluctuations in mRNA pools provide a highly sensitive bioassay for environmental signals that are relevant to microbes
4. Replicated, manipulated experiments fully leverage the value of metatranscriptoimcs for revealing the microbes that percieve a specific environmental change and the metabolic pathways they invoke to respond to it

---

## Metatranscriptomics Workflows:

Cohen NR, Alexander H, Krinos AI, Hu SK, Lampe RH (2024) Marine microeukaryote metatranscriptomics: sample processing and bioinformatic workflow recommendations for ecological applications. _Frontiers in Marine Science_, 9:867007. https://doi.org/10.3389/fmars.2022.867007

Aplakidou E, Vergoulidis N, Chasapi M, Venetsianou NK, Kokoli M, Panagiotopoulou E, Iliopoulos I, Karatzas E, Pafilis E, Georgakopoulos-Soares I, Kyrpides NC (2024) Visualizing metagenomic and metatranscriptomic data: A comprehensive review. _Computational and Structural Biotechnology Journal_, 23:2011-2033. https://doi.org/10.1016/j.csbj.2024.04.060

#### Typical Steps in a Metatranscriptomics Study:

1. Sample collection and RNA extraction
2. cDNA Synthesis
3. Seqeuncing Library Preparation
4. Seqeuncing
5. Data Preprocessing + Quality Control (Trimmomatic, FastQC, MultiQC)
6. Removal of Contaminants and/or Spiked Standards (e.g. non-poly-A RNA and oraganelle mRNA are removed in many marine microeukaryote studies)
7. Assembly
8. Read Mapping
9. Quantification of Trancript Abundance (e.g. when using mRNA standards)
10. Protein Predictions and Orthologous Gene Clustering
11. Taxonomic and Functional Annotations
12. Normalization and Differential Abundance Analysis (e.g. EdgeR, DESeq2)

#### Unique aspects of metatranscriptomics studies:

* The use of mRNA standards to quantify absolute transcript abundance
* It is common to map transcripts directly to genomes you are interested in (e.g. using the MMETSP database)
* You **must** decide what RNA pool to target, and how to treat your sample during library prep

---

How mRNA standards work (figure 1 from Moran et al. 2013): 

<img width="666" alt="Screenshot 2025-04-02 at 6 14 53 PM" src="https://github.com/user-attachments/assets/a4b3fd6d-12a0-4c0b-949c-e40e9c260d2c" />

---

Choosing what RNA pool to target:
1. **ribosomal RNA (rRNA) depletion** - in which rRNA is removed and non-coding RNA and mRNA remain
2. **Polyadenylated (poly-A) RNA selection** - in which mRNA containing poly-A tails (characteristic of eukaryotic mRNA) is selected along with other poly-A-containing non-coding RNA
3. **Sequencing the entire RNA pool** - rRNA will dominate, but this is useful for taxonomic profiling of active microbes (for example, in parallel to eDNA metabarcoding which may recover dormant species or relic/ancient/extracellular DNA)

These two methods differ in their resulting RNA pools, coverage, and quantitaive accuracy, with poly-A selection yielding more protein-coding sequences at a given sequencing depth. 

* rRNA depletion will select for mRNA from **both prokaryotes and microbial eukaryotes**, while poly-A selection will be **biased towards eukaryotes** containing poly-A mRNA 
* Additionally, rRNA depletion better recovers important eukaryotic organelle transcripts without selection bias (plastid expression, e.g. Rubisco, cytochrome oxidase)
* rRNA depletion will will most likely show low recovery of microeukaryotes, and prokaryotes are likely to dominate the majority of the sequencing library

---

Some common metatranscriptomics workflows driven by biological questions (Supp. Figure 1 from Cohen et al. 2021):

<img width="627" alt="Screenshot 2025-04-02 at 6 57 14 PM" src="https://github.com/user-attachments/assets/5f2c56f4-5982-4ea7-96a0-ecb2d43a259e" />

---

Databases and Repositories for Metatranscriptomics Studies (Aplakidou et al. 2024) - resources below we find especially useful are IMG/M, GOLD, KBase, and the NMDC: 

<img width="701" alt="Screenshot 2025-04-02 at 6 35 29 PM" src="https://github.com/user-attachments/assets/2cb0abde-7e27-40e3-856a-f25ed12f1d14" />

---

Some useful visualization concepts and approaches (Figure 2 from Aplakidou et al. 2024): 

<img width="684" alt="Screenshot 2025-04-02 at 6 40 01 PM" src="https://github.com/user-attachments/assets/1107e748-6df0-4c3f-91b1-6a065f3cd466" />

---

Table 3 from from Aplakidou et al. 2024 also has a partial list of tools and software that can be used in -Omics studies. We haven't talked about Geone Browsers much, but you may find JBrowse and Genious especially useful for exploring Genomic and Muti-Omics datasets: 

<img width="671" alt="Screenshot 2025-04-02 at 6 42 48 PM" src="https://github.com/user-attachments/assets/4d5befd5-3c7e-454e-bc7a-2ffb5912030b" />



