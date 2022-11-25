
Some example queries have been given in notebooks [here](https://github.com/ElucidataInc/polly-python/tree/main/Discover)

#### Querying datasets in GEO OmixAtlas
1. To identify datasets belonging to the tissue Breast, disease Breast Neoplasms and organism Homo sapiens

    ```
        query = """SELECT * FROM geo.datasets
                            WHERE CONTAINS(curated_disease, 'Breast Neoplasms')
                            AND CONTAINS(curated_tissue, 'Breast')
                            AND CONTAINS(curated_organism, 'Homo sapiens')
                            """
    ```

2. Fetch datasets from Depmap which has gene dependency information according to CRISPR screening experiments

    ```
        query = """ SELECT dataset_id FROM depmap.datasets
        WHERE
            data_type = 'Gene dependency' AND
            data_repository_version = '2021_Q4' AND
            platform = 'CRISPR Screening'"""
    ```


3. Identify all transcriptome datasets in Hepatocellular Carcinoma disease in Human and Mouse

    ```
        query = """SELECT * FROM geo.datasets
                            WHERE CONTAINS(curated_disease, 'Carcinoma, Hepatocellular')
                            AND (CONTAINS(curated_organism, 'Homo sapiens' OR CONTAINS(curated_organism, 'Mus musculus')
                            AND data_type LIKE '%Transcriptomics%')
                      """
    ```

#### Querying samples in GEO OmixAtlas
1. Get the name of samples, dataset ID and extract_protocol_ch1 where spectroscopy is mentioned in the extract_protocol_ch1

    ```
        query = """SELECT name, src_dataset_id, extract_protocol_ch1 FROM geo.samples
         WHERE LOWER(extract_protocol_ch1) LIKE '%spectroscopy%'"""
    ```

2. IGet the name of disease and number of samples where information about particular disease is curated

    ```
        query = """SELECT curated_disease, COUNT(*) AS count FROM geo.samples 
        GROUP BY curated_disease ORDER BY count DESC """
    ```

#### Querying data matrix in GEO OmixAtlas
1. Fetch data matrix for selected genes for a dataset ID of interest

    ```
        gene = ('hbb-y', 'fth1', 'bbip1', 'actb')
        query = f"SELECT * FROM data_matrices.geo__GSE4230_GPL1261 WHERE 
        LOWER(rid) IN {gene}"
    ```

#### Other Query Examples
1. Select a few feature level metadata for selected genes from Mutation datasets of TCGA where dataset_id contains BRCA

    ```
    query = """SELECT src_dataset_id, disease, protein_position, amino_acids, sequencer, impact, variant_class, consequence, name
        FROM tcga.features AS features
        JOIN (
        SELECT dataset_id AS dataset_id, curated_disease AS disease FROM tcga.datasets WHERE data_type LIKE 'Mutation') AS datasets
        ON features.src_dataset_id = datasets.dataset_id
        WHERE hugo_symbol IN ('TP53','PIK3CA','CDH1','GATA3') AND features.src_dataset_id LIKE '%BRCA%'
        ORDER BY features.src_dataset_id"""
    ```
2. Some cross OmixAtlas querying for Mutation datasets can be found in this [notebook](https://github.com/ElucidataInc/polly-python/blob/main/Discover/Mutation_data_querying%20.ipynb). Here, we have shown query across TCGA, cBioportal and CPTAC OmixAtlas on Polly.