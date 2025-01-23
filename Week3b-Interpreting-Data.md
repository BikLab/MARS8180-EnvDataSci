### Week 3b Bioinformatics Decision Points

_No amount of high-end bioinformatics can compensate for poorly prepared samples and it is therefore imperative that careful attention is given to sample preparation and library generation within workflows, especially those involving multiple PCR steps._

* Quote from Murray et al 2015, From Benchtop to Desktop: Important Considerations when Designing Amplicon Sequencing Workflows - https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0124671

#### 11:10 - 11:30am Bioinformatics Workflows in -Omics studies

Graphics taken from Ben Callahan's lectures at the STAMPS 2024 course at MBL: https://github.com/mblstamps/stamps2024/wiki#17 

<img width="1354" alt="Screenshot 2025-01-22 at 7 10 46 PM" src="https://github.com/user-attachments/assets/8f03a7e2-659f-4c78-b048-a6762d43f25f" />

<img width="1383" alt="Screenshot 2025-01-22 at 7 11 51 PM" src="https://github.com/user-attachments/assets/f2611700-52ad-48f8-88f0-8a129c9fea8b" />

<img width="1399" alt="Screenshot 2025-01-22 at 7 12 15 PM" src="https://github.com/user-attachments/assets/04b21029-122d-4fea-8223-6ed0f9c384e5" />

<img width="1376" alt="Screenshot 2025-01-22 at 7 12 58 PM" src="https://github.com/user-attachments/assets/b483fd49-b02b-42b9-98a4-27be923678b8" />

<img width="1362" alt="Screenshot 2025-01-22 at 7 13 30 PM" src="https://github.com/user-attachments/assets/d37c6ce3-5929-45a0-ab59-b640c1c51012" />

<img width="1355" alt="Screenshot 2025-01-22 at 7 21 53 PM" src="https://github.com/user-attachments/assets/675c8f94-8139-4bc9-8e99-0c770a302730" />

<img width="1171" alt="Screenshot 2025-01-22 at 7 22 55 PM" src="https://github.com/user-attachments/assets/ef76d1be-e762-4769-bc39-a367a55a2c3d" />

<img width="1161" alt="Screenshot 2025-01-22 at 7 23 21 PM" src="https://github.com/user-attachments/assets/557c2974-8706-4bab-af8b-ced42250a4b6" />

<img width="1163" alt="Screenshot 2025-01-22 at 7 24 22 PM" src="https://github.com/user-attachments/assets/a1099931-e2a1-4a86-b516-6978031179fd" />

<img width="1168" alt="Screenshot 2025-01-22 at 7 25 13 PM" src="https://github.com/user-attachments/assets/8d68e10e-3064-4964-83c5-a798bc57f8fa" />

<img width="1155" alt="Screenshot 2025-01-22 at 7 25 42 PM" src="https://github.com/user-attachments/assets/48989ecb-8df0-45ad-a591-396539f6bbd0" />

<img width="1167" alt="Screenshot 2025-01-22 at 7 26 01 PM" src="https://github.com/user-attachments/assets/34846c65-3bfa-4323-85a3-235b07a758c9" />

<img width="1150" alt="Screenshot 2025-01-22 at 7 26 29 PM" src="https://github.com/user-attachments/assets/838f82e3-8282-4167-9f86-267d62c6835c" />

<img width="1138" alt="Screenshot 2025-01-22 at 7 26 54 PM" src="https://github.com/user-attachments/assets/08934931-64bd-4cc7-ade6-d670525d4b40" />

#### 11:30-11:40 Decision points in typical -Omics workflows:

* **Initial data QA/QC** - interpretation of phred scores, choosing quality cutoff, merging criteria (for Paired-End Illumina reads)
* **Clustering or assembly of raw reads** - what algorithm/pipeline to use? Using custom or default software parameters? 
* **QA/QC of clusters/contigs** - Discard singletons? Discard clusters/contigs not meeting certain criteria (e.g. exclude all Amplicon Seqeunce Variants with <50 reads)? Assessment of contamination and/or filtering out potential contaminant seqeunces?
* **Assigning taxonomy or gene function** - What database to use? How to match sequences with taxonomy/functional information (sequenced-based, ontology/hierarcy based, kmer or other method, etc.)?
* **How to filter, categorize, and interpret taxonomy or functional information** - Data mining for a specific species/gene/pathway? Broad exploration of data? Grouping into biologically meaningful categories? Employ statistical analyses, e.g differential abundance across categories/genes/taxa?
* **How to visualize patterns in your data** - SO MANY DECISIONS, this is very hard. Use existing software packages with standard visualizations? More freeform exploration of data using Base R or your own scripts?
* **How to create a scientific narrative through your choices + figures** - this is the culmination of all the above decisions. Are you looking for a specific story (and taking a narrow path of analysis), or do you not know in advance what you should be looking for in your data (and analyzing/visualizing patterns as broadly as possible)? Probably a combination of both, in reality - looking for one story, but keeping your eye out for other patterns.

#### 11:45-11:50am Group Brainstorming

Take 5 minutes and silently brainstorm your most pressing questions on "Bioinformatics Decisions Points" - what things are you struggling with in your own analyses, or what is one area where you need to learn more about to succeed in your own research? 

GDoc Link: INSERT HERE

#### 11:50 - 12:00pm Compiling your study metadata!

Get into pairs, and take 10 minutes to brainstorm a giant list of study metadata that will be relevant to your analysis - this should include things you already have in hand, and things you may need to get from other people (or are results you are waiting on from a lab analysis). This can be anything related to the environment, sample site, time / date / location of sample collection, contextual information, etc.

GDoc Link: INSERT HERE

#### 12:00 - 12:10pm Group reports & discussion

#### 12:10 - 12:25pm Introducing our class metadata file 

Take a peek into "metabarcoding-dataset" mapping and contextual info on GitHub

NOAA Omics data management guide: https://noaa-omics-dmg.readthedocs.io/en/latest/metadata-guidelines.html

Keemei: cloud-based validation of tabular bioinformatics file formats in Google Sheets. Rideout JR, Chase JH, Bolyen E, Ackermann G, González A, Knight R, Caporaso JG. GigaScience. 2016;5:27. http://dx.doi.org/10.1186/s13742-016-0133-6

Keemi Website: https://keemei.qiime2.org/

Keemi Demo metadata file: https://docs.google.com/spreadsheets/d/1_gE_jQcoYGld9aW_dTyE86zdmg1CkNIPHvVJ6CkYvKY/edit?gid=1402180572#gid=1402180572

