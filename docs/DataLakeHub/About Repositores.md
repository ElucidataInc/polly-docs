# 1.3 Index of repositories on Polly Discover

## 1.3.1 GEO

**Introduction**

GEO (Gene Expression Omnibus), managed by the National Center for Biotechnology Information (NCBI) is a public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community.

The three central data entities of GEO are Platforms, Samples, and Series.

Platform - A Platform is composed of a summary description of the array or sequencer. Each Platform record is assigned a unique GEO Identifier (GPLxxx).

1. Sample - A Sample record describes the conditions under which an individual Sample was handled, the manipulations it underwent, and the expression of each element derived from it. Each Sample record is assigned a unique GEO Identifier (GSMxxx). A Sample entity must reference only one Platform and may be included in multiple Series.

2. Series - A Series record links together a group of related Samples and provides a focal point and description of the whole study. Series records may also contain tables describing extracted data, summary conclusions, or analyses. Each Series record is assigned a unique GEO Identifier (GSExxx).

3.  Dataset represents a curated collection of biologically and statistically comparable GEO Samples referring to the same Platform. To enforce this dataset definition, on Polly, each dataset is denoted using a unique ID of the format GSExxxx_GPLxxxx.

For example, GSE100003_GPL15207 would translate to a collection of samples from the series GSE100003 sequenced using the platform GPL15207.

Samples sequenced using X number of Platforms in a particular Series would result in X number of Datasets.

On Polly, you can find all the information corresponding to a dataset in one place (i.e. one GCT file) which includes, study metadata (title, description, author, disease, organism, etc.), sample metadata, and expression data. This saves a considerable amount of time and effort in finding relevant metadata and mapping it to the expression data which can be better spent on the analysis of data.

Moreover, each dataset has been annotated with study metadata fields such as disease, organism, drug, tissue, and dataset ID that can be used to identify relevant dataset(s) on Polly.

**Types of Omics Datasets**

- Bulk Transcriptomics
  -  Microarray
  -  RNA Sequencing


**Usage**

This is the largest curated transcriptomics repository with over 50k datasets. A user can easily identify dataset(s) of relevance and perform analysis. Analysis of Transcriptomics data generally includes comparing specific pairs of samples. The differences may be due to different phenotypes (samples from diseased or healthy tissue, samples with different treatments, samples at different time points undergoing the same treatment, etc.). More commonly, healthy and disease sample groups are compared to discover quantitative changes in expression levels of genes between the two groups, in turn identifying differentially expressed Genes.

Furthermore, the following analyses can be done on top of the differentially expressed Genes:

1. Gene Ontology (GO) and Pathway enrichment analysis
2. Identification of alternative splicing events and Single Nucleotide Polymorphisms (SNPs)
3. Analysis of the protein-protein interaction network

**Level of curation**

Since most of the RNAseq datasets in GEO do not have processed counts and only have raw data in the form of FASTQ files, we did raw data processing using our in-house star alignment pipeline to generate counts. This counts data is then VST normalized using DESEQ package before being wrapped along with the sample metadata as a Dataset.

**Dataset Level** 

We have mapped the following study metadata fields to an ontology so that they remain consistent throughout the repository and querying based on these fields yield appropriate results:
1. Disease
2. Tissue
3. Cell type
4. Cell line
5. Drug
6. Organism

**Sample Level** 

We have mapped the above-mentioned fields in the sample metadata to an onology so that they remain consistent throughout the repository and querying based on these fields yield appropriate results.

We have also deployed our proprietary Machine Learning Model that accurately identifies the samples as Perturbation and Control, which can allow the user to seamlessly analyze a large number of datasets in batches.

**Source**

https://www.ncbi.nlm.nih.gov/geo/


## 1.3.2 LINCS (The Library of Integrated Network-Based Cellular Signatures)
**Introduction**

LINCS program is an initiative by NIH to create a network-based understanding of biology by cataloging the gene expression as well as other cellular processes. When we expose cells to a variety of perturbating agents then it causes a change in gene expression as well as other cellular processes. Developing the network-based approach, it will enable a new understanding of health and disease through an integrative approach the will help to identify the patterns of common network and cellular responses across different types of tissues and cell in response to a broad range of perturbations.

**Types of Omics Datasets**

- Microarray

**Usage**

By generating and making public data that indicates how cells respond to various genetic and environmental stressors, the LINCS project will help us gain a more detailed understanding of cell pathways and aid efforts to develop therapies that might restore perturbed pathways and networks to their normal states.Polly enables the user to query metadata search  across all annotations associated with perturbations, model systems, and signatures.

LINCS database played a crucial role to investigate the reproducibility of the prototypical perturbational assay: quantifying the responsiveness of cultured cells to anti-cancer drugs and influential in the requirement of FAIR data. With additional curations available as part of Polly, facilitates cross comparison in turn ensuring reproducibility.

Identifying the transcription factors (TFs) responsible for observed changes in gene expression is an important step in understanding gene regulatory networks. Enrichment results from these distinct sources are integrated to generate a composite rank that improves the prediction of the correct upstream TF compared to ranks produced by individual libraries.

**Level of curation**

Changes in each cell line measured against treating it with different perturbations which can be drugs or genetic perturbations (CRISPR knockdown perturbations). Each dataset contains the gene expression values for perturbation as well as control vehicle, control untreated and control vector for the respective cell line.  Moreover, these datasets are curated at their dataset and sample level to make them standardized and consistent.

**Dataset Level**

At the dataset level, we tag the metadata of the dataset with ontologies. The fields in the dataset level metadata are
1. disease
2. tissue
3. organism
4. drug
5. author of the study

**Sample Level** 

We have mapped the above-mentioned fields in the samples of the datasets to ontologies so that they remain consistent throughout the repository and querying based on these fields yields appropriate results.

**Source**

Level 3 LINCS data for a gene and drug perturbations has been taken from GEO, using GSE70138 and GSE92742.

* https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138 
* https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
* https://lincsproject.org/LINCS/

## 1.3.3 DEPMAP
**Introduction**

With the advent of knowledge in the molecular basis of cancer, the next question that faced researchers was to establish cause and effect relationship. DepMap consortium is a step in that direction to make the dependency map between the studied alterations in cancer, available drug molecules and physiological processes while enables identifying small molecule sensitivities and predictive biomarkers.

**Types of Omics Datasets**

- Gene Dependency
- RNAi
- Drug Screens
  
**Usage**

After the phenomenal success of the TCGA portal in mankind's fight against cancer, DepMap is to be the most important next step to understand how the molecules (drugs) affect the overall physiological workflow.

DepMap scientists are profiling hundreds of cancer cell line models for genomic information and sensitivity to genetic and small molecule perturbations. The information available as part of results in pooled-cell line chemical-perturbation viability screens has the potential to replace the initial part on drug lead screening. The data for about 4518 compounds can very efficiently be used for novel drug discovery pipelines as well as for approved drug repositioning or repurposing.

CRISPR-Cas9 viability screens are increasingly performed at a genome-wide scale across large panels of cell lines to identify new therapeutic targets for precision cancer therapy. The integration of these data with the already available pool of molecular understanding of cancer recapitulate findings from the individual datasets, provide greater statistical power to cancer- and subtype-specific analyses, unveil additional biomarkers of gene dependency, and improve the detection of common essential genes.

The DepMap dataset is extensively used to perform and predict the interaction between genes. Few therapies target the loss of tumor suppressor genes in cancer. CRISPR-SpCas9 and RNA-interference loss-of-function screens enable the identification of new therapeutic targets associated with genomic loss of tumor suppressor genes.

Additionally, computational biologists, are efficiently using DepMap for developing in-silico methods to detail the in-vivo process. Cell lines are key tools for preclinical cancer research, but it remains unclear how well they represent patient tumor samples. Direct comparisons of tumor and cell line transcriptional profiles are complicated by several factors, including the variable presence of normal cells in tumor samples. Computational tools developed using DepMap data could be used to guide the selection of cell lines that more closely resemble patient tumors and improve the clinical translation of insights gained from cell lines.

**Level of curation**

DepMap repo contains a wide variety of dataset including Gene dependency, Drug screens and RNAi. Though the data is quite structured as far as the DepMap repository is concerned the metadata of cell lines and drugs is available as separate objects. While curating the data, merging of the different information from cell lines, Drugs, probe performance was merged using precise unique IDs to make the holistic information available for the user.

**Dataset Level** 

At the dataset level, we tag the metadata of the dataset with ontologies. The fields in the dataset level metadata are : 
1. disease
2. tissue
3. organism
4. drug
5. description
6. author of the study

**Sample Level**

The data includes information for failed screens which was cleaned as part of curation process.
In order to facilitate quick learning about the cell lines used, Cell line info is made available as column metadata.
For RNAi, we do not currently have a globally standardized reference database, hence the probe Id along with gene symbols were used to annotate the metadata.

**Source**

https://depmap.org/portal/
 


## 1.3.4 Metabolomics
**Introduction**

Metabolomics repository contains all the publicly available metabolomics data sourced from two public repositories - Metabolomics Workbench and Metabolights. Metabolomics Workbench serves as a national and international repository for metabolomics data and metadata. MetaboLights is also a database for Metabolomics experiments and derived information. The database is cross-species, cross-technique. We curate the data from these repositories and make it available in a FAIR manner.

***Datasets in Metabolomics Workbench (MeWork)*** 

Each study on MeWork can have multiple datasets designated by analysis ids. The metadata template of MeWork consists of the following key sections to cover metadata reporting standards recommended by MSI: project, study, experimental design, subjects, treatment, collection, sample preparation, chromatography, analysis and MS/NMR. Each metadata section in turn contains a set of required and optional data fields pertaining to appropriate details about the experiment. For example, the analytical metadata sections include details regarding sample storage conditions, sample preparation and extraction protocols, sample procurement and analytical methods. The dataset id on Polly is framed as STUDYID_ANALYSISID of the corresponding study on MeWork.


***Datasets in MetaboLights***

The two data entities in each MetaboLights Study are MAF and Sample Metadata

Metabolite Assignment File (MAF) - A TSV file containing information about the metabolites investigated in the study. Information regarding database accession IDs , wherein the spectra the metabolite is found and data pertaining to its abundance within the study samples is present in this file. The file name is of the format : m_MTBLSxxx_POS/NEG_Sample_MS/NMR_maf.tsv. Example : m_MTBLS1080_POS_LC-LTQ-MS_metabolite_profiling_v2_maf.tsv

Sample Metadata - The sample information provides all relevant facts about each sample and any controls/standards included in the study. There is ONE sample file corresponding to each MAF file. Sample metadata includes a unique sample name, organism, organism part (for controls use eg. experimental blank and solvent) and sample type (ie. control, QC, experimental sample). Further sample descriptors may include other columns such as Gender, Age, Treatment, etc.  The file name is of the format : a_MTBLSxxx_POS/NEG_Sample_MS/NMR.txt.Example : a_MTBLS1080_POS_LC-LTQ-MS_metabolite_profiling.txt

Study - A study can have multiple MAF and Sample Metadata Files. Each study record is assigned a unique MetaboLights Identifier (MTBLSxxx).A Dataset represents a curated collection of MAF and its corresponding Sample Metadata. To enforce this dataset definition, on Polly, each dataset is denoted using a unique ID of the format MTBLSxxx_m_MTBLSxxx_POS/NEG_Sample_MS/NMR.For example, MTBLS100_m_dwhsaliva_metabolite_profiling_NMR_spectroscopy would translate to a collection of saliva samples from the study MTBLS100 processed using the NMR Spectroscopy protocol. On Polly, you can find all the information corresponding to a dataset in one place (i.e. one GCT file) which includes, study metadata (title, description, author, disease, organism, etc.), sample metadata, and abundance data. This saves a considerable amount of time and effort in finding relevant metadata and mapping it to the abundance data which can be better spent on the analysis of data. Moreover, each dataset has been annotated with study metadata fields such as disease, organism, Ion Mode, tissue, and dataset ID that can be used to identify relevant dataset(s) on Polly.

**Types of Omics Datasets**

- Metabolomics
- Lipidomics

**Usage**

Default application on Polly allows you to perform downstream analysis on single-mode (either positive or negative mode) as well as dual-mode (both positive and negative mode) targeted, semi-targeted (without retention time) and untargeted unlabeled metabolomics data along with insightful visualizations. The app provides a variety of normalization methods, scaling options and data visualization functionalities, thereby allowing an efficient analysis of the data to get actionable insights.

- The application supports data with a simple matrix having samples in the columns and metabolites in the rows.
- It provides different normalization and scaling methods to perform on the data.
- Performs quality checks for internal standards, metabolites, and samples.
- Performs statistical analysis using limma and provides interactive visualizations.
- Performs pathway enrichment analysis and provides pathway visualizations.
- Provides heatmap visualization along with different algorithms like hierarchical clustering, k-means, correlation etc.
- Performs comparative analysis for the different cohort comparisons.

Users can also start a notebook and analyze the datasets using many popular packages like metaboanalyst.

**Level of curation**

Each analysis contains the intensity of metabolites and related metadata for the analysis. Thus the curation for these datasets includes curation at the dataset level.

**Dataset Level** 

At the dataset level, we tag the metadata of the dataset with ontologies. The fields in the dataset level metadata are disease, tissue, organism and drug as well as the description and author of the study.

**Sample Level**

The sample level information in this repository is yet not curated, we are doing both automatic and manual curation for this repository.

**Caveats** 

We discourage comparing multiple datasets in this repository since the experimental conditions  metabolomics can be vastly different.

**Source**

This repository consists of datasets from mainly two sources, Metabolomics Workbench and MetaboLights.

- https://www.metabolomicsworkbench.org/about/index.php 
- https://www.ebi.ac.uk/metabolights/


## 1.3.5 Single Cell
**Introduction**

Single cell sequencing technologies have grown exponentially in the last few years in terms of the number of cells that can be sequenced. With the advent of de-multiplexing platforms like 10x, it is possible to sequence hundreds of thousands of cells in one go. The single cell repository on Polly contains data from these high-throughput platforms(10x, InDrops) as well as datasets generated using older platforms (such as CELSeq, SmartSeq) where the throughput of cells is relatively lower. This is the largest curated single cell RNASeq data repository in the world with over 1500 datasets. Most of the datasets in this repository are sourced from the Gene Expression Omnibus(GEO). It is updated every day to include the recent single cell RNASeq datasets published on GEO. 
This repository also contains singlecell RNASeq datasets from the Human Cell Atlas.

**Types of Omics Datasets**

- Single cell RNA-seq

**Usage**

The use of Single-cell RNA sequencing (scRNA-seq) in understanding cellular biology and disease mechanisms is unparalleled. It allows us to identify and discover complex and rare cell populations, understand the regulatory relationships between genes at the cellular level, and map the trajectory of cellular differentiation. Understanding a biological system in such a precise manner enables us to catch biological signals otherwise missed at the bulk level and design targeted therapies. 

Different types of analyses can be performed on these single cell datasets to answer a variety of biological questions. Clustering and biomarker analysis- to identify the cell types in the population, trajectory inference- to understand the differentiation paths adopted by different sub-populations,  differential gene expression analysis- to identify biomarkers of cell types and genes that differentiate two cell types are some examples of the vast number of analyses that can be done on single cell RNASeq datasets.

These analyses can be performed using various open source packages such as  [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Seurat](https://satijalab.org/seurat/) and [Monocle](https://github.com/cole-trapnell-lab/monocle3) among others. We store the datasets in this repository in the h5ad file format. The h5ad file format is an HDF5 file format that is widely accepted in the single cell sequencing community. It is designed to store large amounts of data, and allow fast querying of parts of a file without accessing the complete file in memory. The dataset h5ad files can be consumed using Polly notebooks, via popular single cell analysis packages like Scanpy and Seurat, which are pre-installed in Polly Notebook environments. They can also be consumed through a GUI-interface, using the [Cellxgene](https://github.com/ElucidataInc/polly-docs/blob/discover_doc_revamp/docs/Data%20Lake%20Revamp.md#2121-proprietary-applications) app or the [single cell visualisation app](https://github.com/ElucidataInc/polly-docs/blob/discover_doc_revamp/docs/Data%20Lake%20Revamp.md#2121-proprietary-applications) hosted on Polly. 

**Level of curation**

A lot of the datasets in this repository are a collection of samples whose cells have been sequenced together. Thus the curation for these datasets includes curation at the study level (title, description, author, disease, organism, etc.), sample level, as well as cell level.  

**Dataset Level**

At the dataset level, we tag the metadata of the dataset with ontologies. The fields in the dataset level metadata are disease, tissue, organism and drug as well as the description and author of the study.

**Sample Level**

We have mapped the above-mentioned fields in the samples of the datasets to ontologies so that they remain consistent throughout the repository and querying based on these fields yields appropriate results. 

**Cell level**

In some datasets we have also manually curated the cell types present in the dataset. This means, that each cell present in the dataset is tagged with a cell type within the context of the study. 

**Source**

- Gene Expression Omnibus - https://www.ncbi.nlm.nih.gov/geo/
- Human Cell Atlas- https://data.humancellatlas.org/

## 1.3.6 TEDDY
**Introduction**

The TEDDY study - The Environmental Determinants of Diabetes in the Young - is looking for the causes of type 1 diabetes mellitus (T1DM). Research tells us that children who get diabetes have certain kind of genes that make them highly susceptible to getting diabetes. However, not all children who are in the high-risk category get diabetes.  It is believed that something happens to "triggers" a child with these risky genes to actually get diabetes. It is the purpose of this study to try and find out what are the triggers that cause children to get diabetes.

The study encompasses results from ~11,000 patients.

**Types of Omics Datasets**

Lipidomics - Case Study [ST001636](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001636) (Lipidomics)

Metabolomics - Case Study [ST001386](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001386) (Metabolomics)

**Usage**

Given the fact that the data ultimately comprises of Lipidomic and Metabolomic expression intensities (but that of young children), they can be easily combined with any other library of Lipidomic or Metabolomic studies (mostly of Adults) to possibly find if any Perturbation group in other studies have similar Lipidomic or Metabolomic features in children and, hence, possibly warn them or even prevent onset of an illness by taking suitable preventive measures to possibly mitigate or even cure the child of the said perturbation.

**Level of curation**

Whenever an instant (a child) is found to be ‘IA’ positive, control instants are then chosen from amongst the remaining samples such that they have similar features (like diet, viral infections, etc). A more in-depth analysis into how do they do the sampling can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/). These pairs are tied together using a Case Index which is unique to each such pair.

For each Case Index, we can have multiple samples under Control and also under Case.

**Sample Level**

The sample [meta-data](ftp://www.metabolomicsworkbench.org/Studies/ST001386_TEDDY_GCTOF_study_design_data_dictionary.xlsx) for each child has the following attributes which can be set by the instant which is shown to be ‘IA’ positive and it is referred to as a ‘case’:

- 'ia_case' -  a Boolean variable: 1 means it is positive for 'IA'
- 'ia_case_ind_1' - it is a case-control pair ID that helps pair up the case instant (this one) with its control instant. Unique for all pairs.
- 'ia_endptage_1' - this is the age in days when the persistent presence of 'IA' was confirmed (2 successive tests).
- ‘sex’ - Gender

The 'IA''s control counterpart also has a similar definition.

Similarly, whenever an instant comes up which is found to be 'T1D' (Type 1 Diabetes) positive, this 'case's control counterpart is chosen from amongst all those instances which are not 'T1D' positive and have similar features. 

**Dataset Level**

The dataset-level meta-data is inherently the same as the sample-level meta-data as the dataset was sliced according to the unique samples that were present.

In this study, given an instant (sample or child) from a class of Perturbation (either IA positive or T1D positive), we also have a set of Control cases that were similar to the Perturbation when it occurred (features like diet history, illness history, visits, stool sample, etc). Now, having dissented up the data on the basis of sample ID (one sample ID per child in the study), Polly allows us to easily pick one perturbation sample and its corresponding set of control samples by simply filtering according to Which sort of perturbation we want (IA positive or T1D positive), and 

Pick a Case Index number that will pull up one Perturbation Instant and its corresponding set of Control Instances.

**Source**

- Case Study [ST001386](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001386)
- Case Study [ST001636](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001636)


## 1.3.7 TCGA
**Introduction**

The Cancer Genome Atlas (TCGA) is a publicly funded project that aims to catalog and discover major cancer-causing genomic alterations to create a comprehensive “atlas” of cancer genomic profiles.
TCGA molecularly characterized over 20,000 primary cancer and matched normal samples spanning [33 cancer types](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/studied-cancers); generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data.
TCGA is cornerstone in cancer biology. It has enabled a deeper understanding of cancer at molecular levels, helped research in science and technology  that has fastened the pace of drug discovery helped patients at the clinic in real-time. TCGA has been influential in R&D processes of many immuno-therapeutics.

**Types of Omics Datasets**

- Transcriptomics (mRNA,miRNA)
- Proteomics
- Copy Number Alteration
- Single Nucleotide Variations & Small Indels
- Methylation

**Usage**

All the available data can be sliced and diced as per user's exquisite research problem to run downstream analysis including but not limited to pathway enrichment, differential expression, TMB, cross-study summary statistics.
Comprehensive curation done for ontologies related to disease, tissue, source, publication, sample (control/perturbation), etc enables easy integration of data from different assays (expression, mutation, etc.) and different cohorts (BRCA, LUAD, etc) for multi-omics integration analysis for putative biomolecular discovery.
In conclusion, TCGA remains the holy grail for cancer researchers and all the novel findings will facilitate diagnosis, treatment, and cancer prevention.

**Level of curation**

Though meticulous efforts have been done at TCGA to structure and standardize the data, none the less some hurdles still exist in accessing and analyzing the data. Polly aims to reduce and where possible eliminate these barriers by standardizing molecular nomenclature (e.g genes will be represented as HGNC ids, and conversions from other formats like Ensemble will be done whenever necessary) and representing sample IDS across the study/Repository.

**Dataset Level**

Comprehensive curation done for ontologies related to the below
1. Cancer Type/Project
2. Disease Type
3. Disease Stage
4. Data type/Assay type
5. Sample Type (Tumor or normal)
6. Gender
7. Vital Status &/or Cause of death
These enable consistency throughout the repository and querying based on the above fields yield appropriate results

**Sample Level**

Sparsity in clinical data has been a consistent problem to deal, polly enables by cleaning up clinical data and adding additional annotations for pathological status, drugs, response, etc. that are the result of research on primary TCGA data.
The [TCGA barcode](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/) is made up of multiple strings representing various information. Clinical data is represented at patient level, but the various assay data are represented at s sample and/or aliquot level.
For eg. Clinical data will have a unique key to be Patient Barcode (tcga-5l-aat0), whereas the transcriptomics (and other assays too) data will have a unique key at the sample level, whether taken from tumor sample or normal site. A single patient can have multiple tumor sample assessments too, example given for a breast (or lung) cancer patient where a patient has developed tumors in both the breasts (or lungs). Additionally, aliquot Barcode is needed to assess the assay-specific information, for e.g the same sample (tumor biopsy) is used to analyze the DNA and transcriptome (also methylation etc depending on the size of initial biopsy) in such case the Sample barcode for the patient is similar across the assays for that given patient and time point. Considering these representations, it becomes a critical task for the user to merge clinical and assay level data. To ease this for the user, in the curation process, the team has parsed all the relevant barcodes available in TCGA and the mappings have been done accordingly. So that for the user wants the task to merge the clinical data (1 patient 1 record) to transcriptomics data for example ( 1 patient 3 records - normal, tumor left breast, tumor right breast) is already done as part of curation and the data on Polly is ready to use for deeper and insightful analysis.
Additionally, it is a known fact that clinical data in TCGA is a bit sparse, so the team has curated additional papers that provide meaningful information. E.g. though TCGA has the drug data, the associated response is not clearly mentioned, additional annotations (from publications/curation packages) have been done in this regard to provide the user with the information not available within TGGA.

**Source**

https://portal.gdc.cancer.gov/
