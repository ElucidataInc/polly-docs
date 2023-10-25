# Single Cell RNAseq Omixatlas- FAQs

## Introduction

Single-cell RNA sequencing has emerged as the technique of choice for researchers trying to understand the cellular heterogeneity of tissue systems under physiological and pathological conditions. Polly provides structured and curated Single-cell RNAseq data in various counts formats as per the need of the research. The data is machine-actionable and analysis-ready for any downstream needs.

Polly can bring in data from a variety of public sources as per the needs as long as the raw counts matrix and metadata are available. Some of the common public sources for Single-cell RNAseq data are:

1. [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
2. [Expression Atlas](https://www.ebi.ac.uk/gxa/home)
3. [Human cell atlas](https://www.humancellatlas.org/)
4. [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell)
5. [Tabula sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/)
6. [Covid-19 Cell Atlas](https://www.covid19cellatlas.org/)
7. [Array Express](https://www.ebi.ac.uk/biostudies/arrayexpress)
8. Publications

### How is a Single Cell RNA-Seq dataset defined on Polly?

A Single-cell RNAseq dataset on Polly represents a curated collection of biologically and statistically comparable samples. A dataset on Polly is a single entity based on its representation in publication or source.


**Dataset ID Nomenclature**
| Source | Source Format | Polly Format |
| ------ | ------ | ------ |
| GEO | For every study, there is a series ID and a platform ID given on GEO <br /> _Platform identifier_ : A Platform record is composed of a summary description of the array or sequencer. Each Platform record is assigned a unique GEO Identifier **(GPLxxx)**. <br />  _Series_ : A Series record links together a group of related samples and provides a focal point and description of the whole study. Each Series record is assigned a unique GEO Identifier **(GSExxx)**. <br /> &nbsp;&nbsp; - A Series record can be a SubSeries or  SuperSeries. <br /> &nbsp;&nbsp; - SuperSeries is all the experiments for a single paper/study and &nbsp;&nbsp; is divided into SubSeries which are different technologies. | Format: GSExxxxx_GPLXXXXX <br /> <br /> All dataset IDs on Polly are GEO subseries |
| Expression Atlas <br /><br /> Human cell atlas <br /><br /> Array Express <br /> | Studies submitted on these databases/portals have unique IDs associated with them <br /> <br /> Eg. E-MTAB-9841, E-HCAD-56, E-GEOD-175929, E-CURD-102 etc. | Datasets for which data is fetched from these sources follow the same nomenclature as that of the source. <br /> <br /> Eg. E-MTAB-9841, E-HCAD-56 etc |
| Single Cell Portal(SCP) | Studies on SCP have a unique ID associated with them <br /> Eg. SCP2331 | Studies fetched from the Single Cell Portal follow the same nomenclature as that of the source. The unique dataset ID format is SCPxxxxx <br /> <br /> Eg.SCP2331 |
| Tabula sapiens | Studies associated with Tabula Sapiens Consortium - Human cell atlas of nearly 500,000 cells from 24 organs are available as dataset ID - ‘_Tabula Sapiens- organ name_’ <br /> <br /> Eg. Tabula Sapiens - Spleen | All 24 datasets on Polly are available with a unique dataset ID format as:<br /><br /> '_TS_organ name_' where <br /><br /> TS stands for ‘Tabula Sapiens’ <br /> <br /> Eg. TS_Tongue |
| Covid-19 Cell Atlas | Studies related to Covid-19 disease: Healthy Donors or Patient Donors <br /> <br /> Dataset Identifiers are unique for each dataset/study. The dataset identifiers can be <br /> &nbsp; - Tissue and/or author name <br /> &nbsp;&nbsp; Eg. Lung parenchyma Vieira Braga et al. <br /> &nbsp; - Key tissue/cell type associated with the study <br /> &nbsp;&nbsp; Eg. Peripheral Blood Mononuclear Cells (PBMCs) <br /> &nbsp; - Other format<br /> &nbsp;&nbsp; Eg. COV028A2 | - For datasets where the associated publication is given: A unique dataset ID as the PMID of the paper <br/><br /> - In case, a publication reference is unavailable, the dataset ID is the same as that on the source |
| Publication |    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   -    | Studies for which data is taken from publication directly (not from the above sources) are available with a unique dataset ID as the PMID of the paper <br /> <br /> Eg. PMID35664061 | 



**Split Datasets** 

Some studies are a collection of smaller studies about different experimental conditions. In such cases, studies, and datasets are split based on specific defined categories. A combined study dataset is available on request. Datasets are split based on the following categories:

1. Disease (Normal will not be separated as another dataset, only diseases will be split)
2. Organism
3. Tissue
4. Assay
   
For split datasets, the nomenclature of dataset ID is as follows:

Polly dataset ID_Split factor, wherein the Polly dataset ID is the format mentioned above for different sources and split factor is the name of the disease, organism, tissue, or Assay.

Eg. GSE156793_GPL24676_Kidney, GSE156793_GPL24676_Liver


## What are the different options within the Single-cell datasets on Polly
Based on the needs of the research, Polly can provide a Single-cell dataset as any of the following outputs

**1. Raw unfiltered counts** <br />
This format contains one H5AD file containing raw counts for unfiltered cells and genes as the data matrix, curated sample level metadata, and normalized feature names. More details are in the subsequent sections. This format is useful if the research needs raw counts to be processed in a custom manner. <br />

**2. Polly processed counts** <br />
This format contains 2 H5AD files.<br />
a. One H5AD file is the same as mentioned in the Raw unfiltered counts section. <br />
b. The second one stores the following <br />
&nbsp;&nbsp;i. processed counts using a consistent pipeline for filtering cells/genes, normalizing the counts, and annotating cell types using the markers from the associated publication. <br />
&nbsp;&nbsp;ii. Raw counts after filtering out the cells and genes as per consistent Polly’s SC pipeline. <br />
&nbsp;&nbsp;iii. Curated sample-level metadata and normalized feature names. <br />
This format ensures the data is post-processed in a consistent format which makes different datasets comparable. This is particularly useful when similar analyses need to be performed across datasets to get an insight. <br />
More details are in the subsequent sections. <br />

**3. Author processed counts** <br />
This format contains 2 H5AD files <br />
a. One H5AD file is the same as mentioned in the Raw unfiltered counts section (only if raw counts are available at source) <br />
b. The second one contains counts processed using the publication’s parameters for filtering cells/genes, normalizing the counts, and annotating cell types using the markers from the associated publication. This H5AD file also contains curated sample-level metadata and normalized feature names. <br />
This format helps replicate the analysis from the associated publication with minimal effort. <br />
More details are in the subsequent sections.<br />


### Details

**1. Raws unfiltered counts** <br /><br /> 
**1.1. Starting Point** 

Single-cell raw count data is typically available in several common file formats across most public resources. For raw counts, data from the source is fetched in the following formats<br />
- Matrix Market (MTX) Format
- Tabular Formats (e.g., CSV, TSV)
- Hierarchical Data Format (e.g. h5ad, suerat, h5)

**1.2 Processing details**

The pipeline for raw counts includes the following steps: <br />

1. **Fetching Raw Counts**: The starting point for the pipeline is getting the raw counts from the source by <br />
 a. Either downloading data directly from the source  - A text file is created with links for downloading <br />
 b. Downloading manually 
Associated barcodes and gene probes are fetched with the raw counts

2. **Preparing Data files**: A JSON file is created for the pipeline to run.<br />
&nbsp;  - A JSON file with all the required input parameters is created. As per the parameters given in the JSON file, further steps in the pipeline are run <br />
&nbsp;  - In case raw counts are downloaded from the source directly, the file is &nbsp;unzipped and data is arranged in a structured way

3. **Creating h5ad**: With the input as the above JSON file, the pipeline starts to create the initial h5ad.<br />
- Programmatic checks are performed to ensure the data matrix contains raw counts only<br />
- If the matrix is not found to be in the correct format (as per the raw count acceptable values), the pipeline stops<br />
- Pipeline accepts the raw counts as per defined the acceptable values, i.e. integer values and >1000<br />
- The index for cell data is sample: barcode and information with respect to each sample: barcode (cell id = sample: barcode) is added<br />
- Feature IDs are converted to Hugo Symbols (or Ensembl IDs where conversion to Hugo symbols is not available or ambiguous)<br />
- Initial h5ad is created<br />
4. **Metadata Curation**:
- The curation pipeline adds curated metadata to the file - sample/cell level metadata is added to the h5ad file<br />
- The curation pipeline generates dataset metadata<br />

5. **QC metrics**: Scanpy QC metrics are calculated and added to the h5ad file<br />

6. **Final h5ad**: The final h5ad file is saved, consisting of QC metrics, curated metadata, raw counts matrix

**1.3 Curation and Metadata details**

On Polly, every dataset and its corresponding samples are tagged with some metadata labels that provide information on the experimental conditions, associated publication, study description and summary, and other basic metadata fields used for the identification of dataset/samples. <br />

There are three broad categories of metadata available on Polly such as:<br />
&nbsp; a) *Source Metadata*: Metadata fields that are directly available from the Source and provided in a Polly-compatible format.<br />
&nbsp; b) *Polly - Curated*: These are curated using Polly’s proprietary NLP-based Polly-BERT models and reviewed by expert bio-curators to ensure high-quality metadata. Some of the Polly curated metadata fields are harmonized using specific biomedical ontologies and the remaining follow a controlled vocabulary.<br />
&nbsp; c) *Standard identifiers* - Identifier metadata tag for datasets, giving information on the basic attributes of that dataset/study.<br /><br />

Metadata for every dataset is available at 3 levels:<br />
&nbsp; I) Dataset-level metadata - General information about the overall experiment/study, for eg, organism, experiment type, and disease under study.<br />
&nbsp; II) Sample-level metadata - Information about each sample related to Eg. Drug, tissue.<br />
&nbsp; III) Feature-level metadata: Information on the genes/features that are consistent across samples.  <br /> 

Details on each metadata field present at the dataset, sample, and feature level are given below<br />

**I Dataset Level Metadata**

**1.4 Output H5AD format and details**
H5AD file for raw counts data offering will be available to the user with the following information: <br />

1. Raw counts `<adata.X>` slot: The raw expression matrix is available in this slot, providing information on the genes and cells as provided by the author.  The raw counts will available be as: <br />
a. Integer counts (with an exception of some datasets from Expression Atlas which may contain fractional values)<br />
b. Values >1000<br />
c. Values not less than 0<br />
d. Sparse matrix<br />

2. Complete Sample Metadata in `<adata.obs>` slot: All the cell/sample level metadata information is available in this slot as per the metadata table given above <br />
a. Polly-Curated Metadata <br />
b. Scanpy QC metrics - QC metrics include n_genes_by_counts, total_counts, total_counts_mt, and pct_counts_mt.<br />
c. Source Metadata<br />
d. Standard Identifier Metadata<br />

3. Processing Details `<adata.uns>` slot: Unstructured metadata on the processing-related details is available in this slot such as: <br />
a. Tools/packages along with their versions used to convert the matrix from source to H5AD file<br />
b. curation model version (for eg. PollyBert_v1)<br />
c. ontology version<br />
d. raw file formats<br />
e. scanpy, anndata version<br />

4. `<adata.uns>` slot: QC Metrics such as mt, n_cells_by_counts, mean_counts, pct_dropout_by_counts & total_counts

**1.5 Report**
A detailed data report is available for each dataset which can be viewed /downloaded using the given link. The following information corresponding to the dataset will be available in each report: <br />

A) General Information: Dataset Summary<br />
- Organism
- Number of cells
- Number of samples
- Disease
- Dataset Source
- Data split: Yes/No
  - In the case of ‘Yes’, the name of the factor used for splitting is given. Eg. Tissue, disease, etc.
- Size of the data matrix
- Annotated cell types from Source

B ) Experiment design
- SC chemistry
- Design

C) Pre-processing
- Reference genome to which the sequenced reads were mapped
- Reference genome annotation used for transcript quantification
- Software pipeline and version used for read mapping and expression quantification

D) Sample/cell Level Metadata Attributes

E) Basic Data QC

F ) Data QA Analysis: This section provides a table showing the data-related checks performed internally along with the final output <br />
&nbsp;&nbsp;i) Metadata QA - Quality checks performed for metadata at the dataset and sample/ &nbsp;&nbsp;cell-level metadata<br />
&nbsp;&nbsp;ii) Data matrix QA

**2. Polly Processed Data** <br /><br />
**2.1. Starting point**

For processing SC datasets using Polly’s standardized pipeline, the starting point is the h5ad file containing raw unfiltered counts created for the raw counts data. The inputs for the processing pipeline are:<br />

Raw h5ad file containing unfiltered counts in X slot, along with sample/cell and feature level metadata (HUGO gene symbols as feature IDs, obs columns specifying batch/sample, and other cell level metadata)<br />

**2.2. Processing details**

The pipeline includes the following steps:<br />

**1. Preparing files**

- Raw h5ad containing unfiltered counts loaded in X slot, along with sample and feature level metadata<br />
- A JSON file is created with all the required input parameters and method choices for the processing workflow<br />

**2.Processing Workflow:** Standard pipeline workflow consists of the following stepsguided by best practices recommendations: <br />

a. **Data Checks**

Data checks are done such as<br />

- X contains raw counts
- X is a sparse matrix in csc_format
- mt,ribo, and hemo tags  are added to genes
- Copy of unfiltered counts matrix from X to raw slot

 b. **Quality Control**

With the gene expression matrix ready (raw counts matrix), the next step is to filter out low-quality cells and genes before further processing. Filtering of cells is done on a per-sample basis, using adaptive cutoffs per dataset<br />

- **Cell filtering**

Following the best practices guide, we use ratios instead of hard thresholds. This adapts cutoffs according to the distribution in each dataset. Outlier cells crossing at least one of the following thresholds are removed (MAD is median absolute deviation)

- `log1p_total_counts > or < 5 MAD`,
- `log1p_n_genes_by_counts > or < 5 MAD`, 
- `pct_counts_in_top_20_genes > or < 5 MAD`, 
- `pct_counts_Mt > 3 MAD or > 8%`
- Gene filtering: genes expressed (counts > 0) in less than 3 cells are removed

c. **QC metrics calculation**

QC metrics are calculated for ribo & hemo genes as well, but genes are not removed based on these metrics.<br />

d. **Doublet detection**

Detection of probable doublets and filtering out predicted doublets from the counts matrix. This is done using the Scrublet tool, which has a Python implementation. [Scrublet is also being used in the EBI SC Expression Atlas pipeline and has been known to perform better against other methods, particularly in terms of speed and scalability (e.g. this review and study)]

Default parameters for Scrublet:
- Doublet rate: 0.06
- Min counts: 2 
- Min cells : 3 
- Min gene variability %: 85 
- Principal components: 30

e. **QC metric re-computation**

Following filtering and doublet removal, the QC metrics are re-computed (using `scanpy.pp.calculate_qc_metrics` function) on the updated filtered matrix and stored in obs/var:<br />

- obs: `'total_counts', 'n_genes_by_counts'`, `'pct_counts_in_top_20_genes'`, `'pct_counts_mt'`, `'pct_counts_ribo'`, `'pct_counts_hb'`

- var: `'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts`

f. **Normalization** 

Normalization is performed to reduce the technical component of variance in the data, arising from differences in sequencing depth across cells. It is performed based on the following choices:<br />

- **Method**: library size/total counts normalization (implemented by `scanpy.pp.normalize` function)
  - **target_sum**: `None` (this is the default option in scanpy and amounts to setting the library size to the median across all cells; best practices guide discusses the use of alternate fixed thresholds like 104 or 106 can introduce over-dispersion)
  - **log1p**: `True`
  - **Scaling**: `True, with zero_center=False` (this will be performed after the HVG step)<br />

Note: Before normalization, the filtered counts matrix is saved to a separate “counts” layer (`adata.layers["counts"])`, and X is updated with the normalized counts matrix. A copy of the normalized (and optionally log1p) counts matrix will also be saved in a separate “normalized_counts” layer, this will preserve the normalized data prior to scaling and batch correction steps<br />

g. **HVG identification**

Feature selection is needed to reduce the effective dimensionality of the dataset, and retain the most informative genes for downstream steps like PCA and neighbor graph<br />

- **Highly Variable Genes (HVG)** are identified using the `scanpy.pp.highly_variable_genes` function with the following settings:
  - HVG method: `"Seurat"` (default, uses the log-normalized counts matrix in X)
  - n_top_genes (# HVGs identified): `2000`
  - batch_key setting: `"sample"` (this corresponds to batch-aware feature selection providing a lightweight form of batch correction to reduce technical differences between samples, see e.g. best practices guide)
  - subset: `False` - the expression matrix is not subsetted to the HVGs. All genes are retained and a highly_variable column gets added to the .var slot

h. **Batch Effect Correction**

Batch effect is checked in the data using “sample” as the batch variable

- In a manual workflow, batch effects in the dataset will be checked visually on a UMAP or tSNE as well as the quantitative metrics of batch effects in the data. 
  - Following quant metrics to be adopted (these are provided through scib and scib-metrics libraries, which were released along with a recent benchmarking study of single-cell data integration methods):<br />
  i. Adjusted rand index comparing batch labels and leiden clusters (ARI)<br />
  ii. Normalized mutual information between batch labels and leiden clusters (NMI)<br />
  iii. Principal component regression variance using batch variable (pcr_batch)<br />
  iv. Scaled graph integration local inverse Simpson’s index (Graph iLISI)<br />
  v. Acceptance rate from batch k-BET test (kbet_accept_rate)<br />

  - For i-iii above, values closer to 0 indicate good batch mixing; for iv-v, values closer to 1 indicate the absence of significant batch effects in the data

  - Batch correction will be done using Scanorama or scVI, both of which return a corrected matrix as well as a batch-corrected low-d embedding. These methods are shown to perform favorably in a recent comprehensive benchmarking study of SC integration methods. [Based on the quality of batch correction and scalability to large numbers of cells, one of these methods would be adopted for further use in the standard pipeline]

- Scanorama: Batch-adjusted expression matrix is obtained by running `scanorama.correct_scanpy(adatas, return_dimred=True)` (adatas is a list of sample-wise anndata objects)

- scVI: By default, the scVI model works with an unnormalized counts matrix and expects counts in “counts” layer. scVI returns the batch-corrected and normalized expression matrix. This is obtained by running the following step after model fitting: `SCVI.get_normalized_expression(library_size = median_library_size, return_numpy=True)` (where median_library_size is the median of library sizes across all cells in the matrix)

- No. of training epochs for scVI VAE model is decided based on dataset size as suggested in best practices guide: `max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])`

  - The resulting matrix from each batch correction method is added as a new layer (`'<method>_corrected'` key) and does not overwrite the X matrix. If the average batch-mixing score is higher for the post-correction matrix compared with the uncorrected matrix, only then, X is updated with the batch-corrected matrix and used for further downstream analysis. The method used on each dataset is also captured in the uns slot (example: `adata.uns["batch_correction_applied"] = {'method': 'scanorama'})`<br />

 Note: Harmony-corrected embedding (output of scanpy wrapper `scanpy.external.pp.harmony_integrate)` would also provided as a separate matrix in the obsm slot ('X_pca_harmony'). However, by default, Harmony does not correct the full expression matrix and only adjusts the PCA representation to mitigate batch differences. 

**i. HVGs are re-computed using the processed, adjusted matrix**

**j. Dimensionality reduction**

This step is necessary to reduce the dimensionality of the data prior to clustering, and to enable visualization for exploratory analysis. Following data embeddings are provided by running the relevant scanpy functions:<br />
- PCA (n_pcs = 50 top PCs ranked by explained variance): `'X_pca'` in obsm, `'PCs'` in varm slot
- Nearest-neighbor graph in PC space: n_neighbors=50, dimensionality of PC space = 40 or number of top PCs explaining 90% variance, whichever is smaller; `'distances'` and `'connectivities'` in obsp slot
- Uniform Manifold Approximation and Projection (UMAP): `'X_umap'` in obsm slot
- t-Distributed Stochastic Neighbor Embedding (t-SNE): `'X_tsne'` in obsm slot

**k. Cell clustering**: Clustering is performed on the processed and dimensionally reduced data matrix to identify (probable) similar groups of cells, possibly representing distinct cell types/states. <br />

- Leiden graph-based clustering is performed with a fixed resolution (resolution = 0.8)
  - The nearest-neighbor graph is constructed in the previous dimensionality reduction step. 
  - Cluster labels are added to the obs slot in a ”clusters” column

**l. Cell type annotation of clusters**: This is done based on the marker genes and cell types from the associated publication, using an automated cell type assignment method -ScType. 

- ScType implementation - Input
  - A tab-separated file is created with
    - raw/author cell types and corresponding lists of marker genes (one cell type per line, no headers)
    - Specification on whether each marker gene is supposed to be present (over-expressed) or absent (under-expressed) in the corresponding cell type.
  - The processed expression matrix in X slot (only the slice of the data corresponding to the supplied marker genes is used for scoring cell types).

 - ScType implementation - Output
  - ScType returns an enrichment score matrix of size `num cell types x num clusters`, based on which each cluster gets assigned the raw cell type with the highest score. 

- Annotated Tags:
  - Following sample-level metadata columns get added to obs (these would be uniform across all cells in each cluster):<br />
i. polly_curated_cell_type: Raw cell type associated with the cell<br />
ii. polly_curated_cell_ontology:  Cell type associated with the cell, curated with standard ontology<br />
iii. marker_gene_present: Raw cell type and corresponding present (over-expressed) marker genes corresponding to the assigned cell type against each cell<br />
iv. marker_gene_absent: Raw cell type and corresponding absent (under-expressed) marker genes corresponding to the assigned cell type against each cell<br />

   - Following outputs of the cell annotation step are saved to the uns slot of the anndata object:<br />
i. A table of raw/author cell type predictions by cluster, providing the corresponding ScType scores and confidence values, along with the subset of marker genes for the assigned cell type which are differentially expressed in the annotated cluster.<br />
ii. Marker genes dictionary used for the cell type annotation step ({cell_type:[markers]})<br />

In addition to the above outputs, the differentially expressed genes from one-vs-all comparison per cluster, calculated at the above step, would also get automatically saved to the .uns slot (as `'rank_genes_groups'`)<br />

**m. Saving the data**

The parameter dictionary will be saved to the uns slot, after adding a time stamp as the “time_stamp” key to the params dict. The X matrix will be converted to csc_matrix format and the processed anndata object will be saved in H5AD format and pushed to Polly’s omixatlas. Further, details of the processing performed (parameters/method choices, QC metrics, and associated plots) are bundled as a separate comprehensive HTML report and provided with the data file.

 **2.3. Curation and Metadata details**

Metadata curation for Polly annotated data is the same as raw counts data (given above). However, there are some additional fields that are present for every polly-processed dataset in addition to the existing fields available for raw unfiltered data.<br />

Details on each metadata field present at the dataset, sample, and feature level are given below<br />

**I dataset Level Metadata**




