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

All GEO RNA-Seq datasets on Polly are processed using the Kallisto Pipeline. The data is processed using the following reference genome, annotation, and complementary DNA sequence data from Ensembl release 107 for each organism. However, approximately 12% of the datasets have been processed with the Ensemble release V90. These will be reprocessed in future based on the Ensemble release V107.

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

![](RackMultipart20230112-1-ep1uzj_html_e95051302e35596c.png)

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
14. Each sample is then annotated with relevant metadata using our custom curation models for fields like disease, tissue, cell line, drug etc.
15. If requested, the counts matrix is normalized using DESeq2 VST (Variance Stabilizing Transformation).
16. GCT having Raw Counts is pushed to the Omix Atlas - Bulk RNASeq OmixAtlas.

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

On Polly, all datasets and the corresponding samples are available with six standard metadata fields which are curated using Polly's proprietary NLP-based Polly-BERT models and harmonized with specific biomedical ontologies. Such fields are Organism, Disease, Tissue, Cell Line, Cell Type and Drug. In addition to these fields, other metadata fields are also available which are extracted from the source metadata. Metadata is available at 3 levels:

- Dataset-level metadata - General information about the experiment, subject and transformations, for eg, organism, experiment type, and disease under study.
- Sample-level metadata - Captures information for each sample. Eg. Drug, tissue
- Feature-level metadata: Provides molecular/gene information that is consistent across samples.

### How frequently data is added to the Bulk RNASeq OA?

New datasets will be added to the source Bulk RNASeq OA at a frequency of once per week.

### Can raw data be accessed for a sample?

Raw data can be accessed through the SRA links available in columns starting with the name "relation".

### Are source fields retained for a dataset?

All source fields are captured in the GCT file which users can visualize through applications. However, fields available for searching through UI and Polly-python are restricted by schema due to performance constraints.

### Are external data requests directly added to Bulk RNAseq OA?

Users can request Bulk RNASeq datasets from GEO which are not available in the source Bulk RNASeq OmixAtlas to be added to the destination atlas. Upon request, these datasets will be added directly to the user's destination Atlas for which the TAT is 4 days for up to 10 datasets and 10 days for up to 100 datasets. The scope of the data request is as follows:

1. Raw data files (fastq) are available on GEO
2. The dataset belongs to the Bulk RNA-Seq data type (not Single Cell)
3. Standard metadata fields at the dataset and sample level will be auto curated
