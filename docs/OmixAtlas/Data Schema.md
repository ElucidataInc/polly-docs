
## What is the standardized schema of OmixAtlas?

Data schema: the data available within OmixAtlas is curated within defined indexes on the basis of the information it contains. These indexes are:

- Dataset-level metadata (index: files): Contains curated fields like drug, disease, tissue organism, etc., for each dataset.
- Sample-level metadata (index: gct\_metadata, h5ad\_metadata, and biom\_metadata): Contains curated fields like cell lines, experimental design, etc., for each sample.
- Feature level metadata (gct\_row\_metadata, h5ad\_data, and biom\_data): Contains the gene/molecule symbol along with the feature intensity for each sample.
- Variant-related data (index: variant\_data): Contains the schema for variant-related information present in vcf files
