
## 1. What is the standardized schema of OmixAtlas?

Data schema: the data available within OmixAtlas is curated within defined indexes on the basis of the information it contains. These indexes are:

- Dataset-level metadata (index: files): Contains curated fields like drug, disease, tissue organism, etc., for each dataset.
- Sample-level metadata (index: gct\_metadata, h5ad\_metadata, and biom\_metadata): Contains curated fields like cell lines, experimental design, etc., for each sample.
- Feature level metadata (gct\_row\_metadata, h5ad\_data, and biom\_data): Contains the gene/molecule symbol along with the feature intensity for each sample.
- Variant-related data (index: variant\_data): Contains the schema for variant-related information present in vcf files

## 2. List of curated field on Polly:

There are standard metadata fields that are a part of product offering for every dataset available on Polly. These fields follow a particular ontology system as shown below.

| Field name | Ontology | Description |
| --- | --- | --- |
| curated\_disease | [MeSH](https://www.ncbi.nlm.nih.gov/mesh/1000067) | At the sample-level, the disease is specific to the sample under study. At the dataset-level, this field captures all disease mentions from the text. This includes mention of related diseases, diseases derived from cell-line mentions that might or might be directly under study in the experiment, and all diseases from sample-level metadata.|
| curated\_organism | [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) | This field is usually present only at the dataset-level. |
| curated\_tissue | [Brenda Tissue Ontology (BTO)](https://www.ebi.ac.uk/ols/ontologies/bto) | At the sample-level, the tissue is specific to the sample under study. At the dataset-level, this field captures all tissue/organ mentions from the text. This includes mention of related tissues, that might or might be directly under study in the experiment, and all tissues from sample-level metadata.|
| curated\_cell\_line | [Cellosaurus](https://www.cellosaurus.org/) | At the sample-level, the cell-line is specific to the sample under study. At the dataset-level, this field captures all cell-line mentions from the text. This includes mention of related cell-lines, that might or might be directly under study in the experiment, and all cell-lines from sample-level metadata.At the sample-level, the cell-type is specific to the sample/single cell under study. At the dataset-level, this field captures all cell-type mentions from the text. This includes mention of related cell-types, that might or might be directly under study in the experiment, and all cell-types from sample-level metadata|
| curated\_cell\_type | [Cell Ontology (CL)](https://www.ebi.ac.uk/ols/ontologies/cl) | At the sample-level, the cell-type is specific to the sample/single cell under study. At the dataset-level, this field captures all cell-type mentions from the text. This includes mention of related cell-types, that might or might be directly under study in the experiment, and all cell-types from sample-level metadata.|
| curated\_drug | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | At the sample-level, the drug(s) is specific to the sample under study. At the dataset-level, this field captures all drug mentions from the text. This includes mention of related drug(s), that might or might be directly under study in the experiment, and all drugs from sample-level metadata. |

| **Field Name** | **Metadata Level** |
| --- | --- |
| dataset\_id | Dataset |
| description | Dataset |
| data\_type | Dataset |
| curated\_organism | List |
| curated\_disease | Dataset |
| curated\_tissue | Dataset |
| curated\_drug | Dataset |
| curated\_cell\_line | Dataset |
| curated\_cell\_type | Dataset |
| dataset\_source | Dataset |
| publication | Dataset |
| curation\_version | Dataset |
| total\_num\_samples | Dataset (for single cell datasets) |
| sample\_id | Sample |
| curated\_disease | Sample |
| curated\_cell\_line | Sample |
| curated\_tissue | Sample |
| curated\_drug | Sample |
| curated\_cell\_type | Sample |
| curated\_genetic\_mod\_type | Sample |
| curated\_genetic\_modified\_gene | Sample |

The UI filters that are ontology-backed are disease, tissue and cell-line.

| **repo\_name** | **Data Type** |
| --- | --- |
| **hpa** | transcriptomics, gene expression reliability |
| **gdc** | mirna expression, transcriptomics, copy number variation |
| **geo** | transcriptomics, raw counts transcriptomics, snp array |
| **enterprise\_atlas** | transcriptomics, single cell, mutation, mirna |
| **cptac** | proteomics, phosphoproteomics, transcriptomics, copy number variation, mutation, mirna expression, acetylproteomics, methylation, lipidomics, metabolomics |
| **liveromix\_atlas** | transcriptomics, mutation, metabolomics, single cell, proteomics, lipidomics, mirna, drug screens, gene dependency, gene effect |
| **cbioportal** | mutation, copy number variation, fusion, transcriptomics, methylation |
| **depmap** | gene dependency, drug screens, gene effect, rnai |
| **pcd** | drug response |
| **lincs** | transcriptomics |
| **metabolomics** | metabolomics, lipidomics, single cell |
| **immport** | lab measurement, proteomics, titer, pcr, cytometry |
| **ukbiobank** | gwas |
| **sc\_data\_lake** | single cell |
| **tcga** | transcriptomics, mirna, copy number variation, mutation, methylation, proteomics |
| **gtex** | transcriptomics |
| **gnomad** | gwas |
| **rcsb** | structural biology |
| **teddy** | lipdomics, metabolomics |
| **PRIDE** | |
| **Xena** | |
| **HugeAMP** | GWAS |
| **OpenGWAS** | GWAS |
| **ChemOA** | |
