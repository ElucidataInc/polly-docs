# Introduction

Spatially resolved transcriptomics (SRT) is a cutting-edge scientific method that merges the study of gene expression with precise spatial location within a tissue. This revolutionary approach allows researchers to visualize the spatial distribution of RNA transcripts within a tissue sample, essentially mapping out where each gene is being expressed. Traditional transcriptomics analyses gene expression in bulk tissue samples or single cells but lacks spatial information. 

Polly can bring in SRT data from a variety of public sources as per the needs as long as the raw counts matrix, spatial coordinates, imaging data, and metadata are available. Some of the common public sources for Spatial Transcriptomics data are:

1. [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
2. [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell)
3. [Zenodo](https://zenodo.org/records/4739739)
4. [CZI-CellxGene](https://zenodo.org/records/4739739)
5. Publications



# How is a spatial transcriptomics dataset on polly defined?

Spatial Transcriptomics Datasets on Polly represent a curated collection of biologically and statistically comparable samples. Each dataset is denoted using a unique ID, whose nomenclature depends on the public database it was ingested from. 

Spatial datasets are stored in a sample-wise fashion, i.e. multiple samples from a series are hosted as individual datasets. This is done as most spatial studies analyze individual samples independently due to lack of a straightforward way to combine spatial coordinates and imaging data from distinct tissue sections. 

| Source | Source Format | Polly Format |
| ------ | ------ | ------ |
| GEO |For every study, there is a series ID and a sample ID given on GEO. <br> **Series** - A Series record links together a group of related samples and provides a focal point and description of the whole study. Each Series record is assigned a unique GEO Identifier (GSExxx). A Series record can be a SubSeries or  SuperSeries. SuperSeries is all the experiments for a single paper/study and is divided into SubSeries which are different technologies.</br> <br> **Sample** - A Sample record describes the conditions under which an individual Sample was handled, the manipulations it underwent, and the abundance measurement of each element derived from it. Each Sample record is assigned a unique and stable GEO accession number (GSMxxx) </br>| Format: GSExxxxx_GSMXXXXX |
| Single Cell Portal(SCP) | Studies on SCP have a unique ID associated with them.  Eg. SCP2331 | Studies fetched from the Single Cell Portal follow the same nomenclature as that of the source. The unique dataset ID format is SCPxxxxx  Eg.SCP2331 |
| Zenodo | Zenodo hosts “records”, which are the basic entities used to share and preserve a digital research object (datasets, publications, software, poster, presentations etc). Each record on Zenodo has a unique Digital Object Identifier (DOI) associated with them Eg. 4739738 | Studies fetched from Zenodo follow the same nomenclature as that of the source. The unique dataset ID format is for ex. 4739738 |
| CZI-CellxGene | Datasets on CZI have a unique alphanumeric ID associated with them. Eg. 60358420-6055-411d-ba4f-e8ac80682a2e | Studies fetched from CZI follow the same nomenclature as that of the source. The unique dataset ID format is for Eg. 60358420-6055-411d-ba4f-e8ac80682a2e |
| Publication |   -    | Studies for which data is taken from publication directly (not from the above sources) are available with a unique dataset ID as the PMID of the paper. Eg. PMID35664061 | 


# What are the different options within the Spatial datasets on Polly?

Currently, Polly provides Spatial datasets (generated using the 10X Visium technology) as any of the following outputs. To be noted, the sampled entities for Visium datasets are “spots” (which contain 10-30 single cells) rather than individual cells.  

**Raw unfiltered counts**

This format contains one H5AD file containing raw counts for unfiltered spots and genes as the data matrix, spatial coordinates of the spots, H&E image of the tissue section used for generating the data, curated sample level metadata, and normalized feature names. More details are in the subsequent sections. This format is useful if the research needs raw counts to be processed in a custom manner.

**Custom processed counts**

This format contains 2 H5AD files.

a. One H5AD file is the same as mentioned in the Raw unfiltered counts section 

b. The second one stores the following in addition to fields in the raw data

- Processed counts using a pipeline for filtering spots/genes, normalizing the counts must be as per requirement
- Raw counts after filtering out the spots and genes as per requirements.
- Normalized and log1p transformed counts.
- Results from dimensionality reduction and clustering of spots
- Results from spatially variable genes (SVG) analysis
- Results from cell type deconvolution using associated reference single cell data.
- Curated sample-level metadata and normalized feature names.

The pipeline uses default standard parameters at every processing step which can be customized as per requirement. However, methods used for SVG and cell-type deconvolution are currently fixed.


## Details

### 1. Raws unfiltered counts

**1.1. Starting Point** 

To process spatial transcriptomics datasets using Polly's pipeline, the starting point is a counts matrix in h5 format, spatial coordinates (of the spots) provided in csv format, both low-resolution and high-resolution H&E images of the tissue section used for data generation, a scalefactors JSON file for the mapping of spot coordinates to image pixels, curated sample-level metadata, and normalized feature names. 

**1.2 Processing details**

The pipeline for raw counts includes the following steps: 

1. **Fetching Raw Counts**: The starting point for the pipeline is getting the raw counts and spatial data from the source by

   - Either downloading data directly from the source  - A text file is created with links for downloading 
   
   - Downloading manually

2. **Preparing Data files**: A JSON file is created for the pipeline to run.

   - A JSON file containing all necessary input parameters is generated. Subsequent steps in the pipeline are executed based on the parameters specified in this JSON file.

3. **Creating h5ad**: With the input as the above JSON file, the pipeline starts to create the initial h5ad.


   - Programmatic checks are performed to ensure that counts data, spatial coordinates and required image files (tissue_hires and lowres images, scalefactors.json) are available. 
   
   - Check that the data matrix contains raw counts only
   
   - If the matrix is not found to be in the correct format (as per the raw count acceptable values), the pipeline stops
   
   - Pipeline accepts the raw counts as per defined the acceptable values, i.e. integer values 
   
   - The index for spot data is “cell id” which is in the form of “sample id : spot barcode” and information with respect to each “cell id” is added
   
   - Feature IDs are converted to Hugo Symbols (or Ensembl IDs where conversion to Hugo symbols is not available or ambiguous)
   
   - Initial h5ad is created 


4. **Metadata Curation**:


   - The curation pipeline adds curated metadata to the file - sample/cell level metadata is added to the h5ad file
   - The curation pipeline generates dataset metadata


5. **QC metrics**: Scanpy QC metrics are calculated and added to the h5ad file

6. **Spatial Embedding for visualization with cellXgene:** A X_spatial embedding is added to the h5ad obsm slot. This is a copy of the spatial embedding which is added automatically when creating the initial h5ad.

7. **Final h5ad**: The final h5ad file is saved, consisting of QC metrics, curated metadata, raw counts matrix, spatial coordinates and H&E images.


**1.3 Curation and Metadata details**

On Polly, every dataset and its corresponding samples are tagged with some metadata labels that provide information on the experimental conditions, associated publication, study description and summary, and other basic metadata fields used for the identification of dataset/samples.

There are three broad categories of metadata available on Polly such as:

a) *Source Metadata*: Metadata fields that are directly available from the Source and provided in a Polly-compatible format.

b) *Polly - Curated*: These are curated using Polly’s proprietary NLP-based Polly-BERT models and reviewed by expert bio-curators to ensure high-quality metadata. Some of the Polly curated metadata fields are harmonized using specific biomedical ontologies and the remaining follow a controlled vocabulary.

c) *Standard identifiers* - Identifier metadata tag for datasets, giving information on the basic attributes of that dataset/study.

Metadata for every dataset is available at 3 levels:

I) Dataset-level metadata - General information about the overall experiment/study, for eg, organism, experiment type, and disease under study.

II) Sample-level metadata - Information about each sample related to Eg. Drug, tissue.

III) Feature-level metadata: Information on the genes/features that are consistent across samples. 


**1.4 Output H5AD format and details**

H5AD file for raw counts data offering will be available to the user with the following information:

1. Raw counts  `<adata.X>` slot: The raw expression matrix is available in this slot, providing information on the genes and cells as provided by the author.  The raw counts will available be as: 


   - Integer counts (with an exception of some datasets from Expression Atlas which may contain fractional values)
   - Values not less than 0
   - Sparse matrix 

2. Spatial Coordinates data: Columns in_tissue, array_row, array_col are available in the `<adata.obs>` slot. 


   - in_tissue column: indicates whether the spot falls under the tissue section
   - array_row and array_col: provide location of each spot on the Visium array
   - Pixel coordinates for the spots (for plotting spots on the image) are added as a 2D array to the <adata.obsm> slot (key: spatial) 

3. A “X_spatial” key is added to `<adata.obsm>` slot. This contains exactly same info as adata.obsm[“spatial”] which is automatically added upon creating the h5ad. It is required for compatibility with downstream visualization packages.


4. Imaging Data: Image data is added to the `<adata.uns>`. The “spatial” key in the uns slot holds a dictionary containing lowres, hires image arrays and scalefactors used for mapping spot coordinates onto the lowres/ hires images.  


5. Complete Sample Metadata in `<adata.obs>` slot: All the spot/sample level metadata information is available in this slot as per the metadata table given above


   - Polly-Curated Metadata
   - Scanpy QC metrics - QC metrics include n_genes_by_counts, total_counts, total_counts_mt, and pct_counts_mt.
   - Source Metadata
   - Standard Identifier Metadata


6. Processing Details `<adata.uns>` slot: Unstructured metadata on the processing-related details is available in this slot such as:


   - Tools/packages along with their versions used to convert the matrix from source to H5AD file
   - ontology version
   - raw file formats
   - scanpy, anndata version


7. `<adata.uns>` slot: QC Metrics such as mt, n_cells_by_counts, mean_counts, pct_dropout_by_counts & total_counts




