# Bulk RNA sequencing FAQs

RNA sequencing (RNA-seq) has become a powerful method in transcriptomics, allowing highly accurate gene expression quantification. It is providing researchers with visibility into previously undetected changes occurring in disease states, in response to therapeutics, under different environmental conditions, and across a broad range of other study designs.

### What is a Bulk RNASeq dataset on Polly?

A dataset on Polly represents a curated collection of biologically and statistically comparable samples. All the GEO bulk RNA-Seq datasets on Polly are denoted using a unique ID. This dataset ID follows the GEO record identifier format, consisting of a series ID and a platform ID.

- Platform record - A Platform record is composed of a summary description of the array or sequencer. Each Platform record is assigned a unique GEO Identifier ( **GPLxxx** )
- Series record - A Series record links together a group of related samples and provides a focal point and description of the whole study. Each Series record is assigned a unique GEO Identifier ( **GSExxx** ). A Series record can be of two types - SubSeries and SuperSeries. For bigger experiments, there are both SubSeries and SuperSeries. SuperSeries is all the experiments for a single paper/study. SuperSeries can be divided into SubSeries which are different technologies. In the GEO Transcriptomics OmixAtlas, we only have Subseries records.

For example, dataset ID 'GSE189190\_GPL25947\_raw' would translate to:

- sequenced using the **platform ID- GPL25947**
- from the **Series GSE189190**
- Here \_ **raw** signifies 'raw counts' for datasets

### In which form are the files of Bulk RNASeq datasets available on Polly?

All datasets are available with raw counts data matrix and associated metadata. The datasets are available in GCT (Gene Cluster Format) file format which allows for storing sample metadata as well as the expression data in a single file.

### Which pipeline is used for processing the RNASeq data?

All Bulk RNA-Seq Datasets on Polly are processed using a Kallisto pipeline. The data is processed with the following reference genome, annotation, and complementary DNA sequence data from Ensembl release 107 for each organism. 

**Note**: approximately 12% of datasets on the Bulk RNA-Seq OmixAtlas are currently processed with the Ensemble release V90. These will be reprocessed based on the Ensemble release V107 in future versions

1. **Homo Sapiens** Ensembl release 107, 90
  1. Genome sequence (fasta)
  2. Gene annotation set (GTF)
  3. cDNA sequences (fasta)
2. **Mus musculus** Ensembl release 107, 90
  1. Genome sequence (fasta)
  2. Gene annotation set (GTF)
  3. cDNA sequences (fasta)
3. **Rattus norvegicus** Ensembl release 107, 90
  1. Genome sequence (fasta)
  2. Gene annotation set (GTF)
  3. cDNA sequences (fasta)

**Process flow**

![Process flow](../img/OmixAtlas-Images/1_1.png) <center>**Figure 1.** Process Flow</center>

#### Details of the processing steps:

1. Detect organisms and fetch relevant genome, annotation, and complementary DNA sequence data from Ensembl.
2. Download the transcriptome sequencing data (.sra files) from SRA using sratoolkit prefetch / AWS S2 URI if publicly available.
3. Validate the downloaded .sra file using vdb-validate.
4. Identify if the SRA data is (single-end) or (paired-end)using fastq-dump. Both single-end (SE) and paired-end (PE) sequencing data are processed with the exclusion of color-space sequence data.
5. Extract fastq files with parallel-fastq-dump.
6. Perform basic quality control checks on the .fastq reads using FastQC. (Diagnose basespace / colorspace, quality encoding, read length)
7. Trim Bases with phred quality \<10 on the 3′ ends and discarded reads shorter than 18 nucleotides using Skewer.
8. Adapter sequences at the 3′ end are detected using Minion.
9. If the predicted adapter sequence is not present in the genome and exceeds a frequency of 2.5% then the adapter sequences are clipped using Skewer.
10. Adapter contamination detection using bowtie and clipping using a skewer.
11. Transcript-level expression counts are generated using Kallisto by mapping all the reads that pass quality control to the genome. Command: "kallisto quant" . All counts are reported on the gene level by taking a simple sum of Transcript-level counts. (NOTE: Kallisto pseudo counts are rounded to integer values)
12. For every SRR accession, the generated counts are collected into a single (.gct) file and multiple SRR counts per GSM ID (sample) are aggregated.
13. At the feature level, the Ensembl gene IDs are mapped to the respective HGNC symbol, MGI Symbol or RGI symbol. Counts for duplicate genes are dropped using Mean Average Deviation Score.
14. Each sample is then annotated with relevant metadata using Polly’s curation models for a standard fields **disease, tissue, cell line, drug, cell type, organism.** 
15. If requested, the counts matrix is normalized using DESeq2 VST (Variance Stabilizing Transformation).
16. GCT having Raw Counts are pushed to the Omix Atlas - Bulk RNASeq OmixAtlas.

### Tools Used for the processing:

| | **Tool** | **Task** | **Usage** |
| --- | --- | --- | --- |
| 1 | GEOparse | Query GEO and fetch sample IDs (GSMs). | |
| 2 | pySRAdb | Query SRA and fetch run IDs corresponding to sample IDs (GSMs) and create GSM: SRR mappings. | |
| 3 | SRA toolkit | Download SRA files | prefetch SRRXXXXXX |
| 4 | SRA toolkit |Validate downloaded SRA files | vdb-validate |
| 5 | SRA toolkit |diagnose single or paired-end | fastq-dump |
| 6 | SRA toolkit |dump fastq | parallel-fastq-dump, |
| 7 | FastQC | Diagnose basespace / colorspace, quality encoding, read length | fastqc |
| 8 | parallel-fastq-dump | Rapid decompression of sequence data from .sra files | parallel-fastq-dump |
| 9 | Minion | 3' adapter detection | minion search-adapter |
| 10 | Bowtie2 | Adapter contamination detection | bowtie2 |
| 11 | Skewer | <ul><li>- 3' quality trimming </li><li> Adapter clipping </li><li> 5' trimming </li></ul> | skewer |
| 12 | FASTX-Toolkit | Progressive 5' trimming | fastx\_trimmer |
| 13 | Kallisto | Transcript-level mapping | Kallisto quant |
| 14 | custom script (make GCT) | Collect transcript counts, sample metadata and make counts matrix then make a GCT file | |
| 15 | GEOtron | Curate sample and data-set level information and attach it to the GCT file | |

### **Which fields are curated for bulk RNASeq data Data?**

Polly's NLP-based curation models are used to curate all datasets and their corresponding samples, which are then harmonised using specific biomedical ontologies. These fields include - Organism, Disease, Tissue, Cell Line, Cell Type and Drug. 
In addition to these fields, various other metadata fields are captured from the source publication and annotated at 3 levels:

- Dataset-level metadata - General information about the experiment, subject and transformations, for eg, organism, experiment type, and disease under study.
- Sample-level metadata - Captures information for each sample. Eg. Drug, tissue
- Feature-level metadata: Provides molecular/gene information that is consistent across samples.

### How frequently is the Bulk RNASeq Atlas updated?

New datasets will be added to the source Bulk RNASeq OA at a frequency of once per week.

### Can raw data be accessed for a sample on Polly?

Raw data can be accessed through the SRA links available in metadata tables on the OmixAtlas UI. Users can navigate to the column name - "relation" to access the link.

### Are all metadata fields from the source retained?

All metadata fields from the source publication are captured in the GCT file and can be visualized through applications. However, only the fields described in the OmixAtlas schema are available for querying through the UI and Polly-python.

### Can users request for datasets on Polly?

Bulk RNASeq datasets from GEO which are not available on Polly can be requested to be added to a user's OmixAtlas. Upon request, these datasets will be added directly to your atlas within 4 days. The following requests are in scope:

1. Raw data files (fastq) are available on GEO
2. The dataset belongs to the Bulk RNA-Seq data type.
