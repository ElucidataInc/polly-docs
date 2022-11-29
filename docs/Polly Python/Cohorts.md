::: polly.Cohort
    options:
      show_source: false

## Examples

### In TCGA

```py
query = someSQLquery
results=omixatlas.query_metadata(query)
```

    Query execution succeeded (time taken: 2.13 seconds, data scanned: 0.244 MB)
    Fetched 123 rows



```py
dataset_ids = results['dataset_id'].tolist()
cohort1.create_cohort("/import","tcga_data","Proteomics datasets","tcga", dataset_ids)
```

    INFO:root:Cohort Created !


    Initializing process...


    Verifying Data: 100%|██████████| 123/123 [00:11<00:00, 10.71it/s]
    Adding data to cohort: 100%|██████████| 123/123 [00:14<00:00,  8.72it/s]
    Adding metadata to cohort: 100%|██████████| 123/123 [00:11<00:00, 10.25it/s]
    INFO:root:'123' dataset/s added to Cohort!


```py
dataset_metadata = cohort1.merge_data("dataset")
display(dataset_metadata.head())
```

```py
All_Metadata_col = cohort1.merge_data("sample")
print("\nColumns/Datasets information")
display(All_Metadata_col.head())
```

```py
df_real = cohort1.merge_data("data_matrix")
print("\nData Matrix")
display(df_real.head())
```
### In GEO

```py
dataset_ids = results['dataset_id'].tolist()
cohort1.create_cohort("/import","geo_data","Transcriptomics datasets","geo", dataset_ids[0])

for i in dataset_ids[1:]:
    cohort1.add_to_cohort("geo", i)
```

    INFO:root:Cohort Created !


    Initializing process...
    Adding data to cohort...
    Adding metadata to cohort...


    INFO:cmap_logger:Reading GCT: /import/geo_data.pco/geo_GSE120746_GPL18573.gct
    INFO:root:'18' sample/s added to Cohort!


    Initializing process...
    Adding data to cohort...
    Adding metadata to cohort...


    INFO:cmap_logger:Reading GCT: /import/geo_data.pco/geo_GSE62642_GPL16791.gct
    INFO:root:'14' sample/s added to Cohort!


    Initializing process...
    Adding data to cohort...
    Adding metadata to cohort...


    INFO:cmap_logger:Reading GCT: /import/geo_data.pco/geo_GSE68719_GPL11154.gct
    INFO:root:'73' sample/s added to Cohort!


```py
dataset_metadata = cohort1.merge_data("dataset")
display(dataset_metadata.head())
```

```py
All_Metadata_col = cohort1.merge_data("sample")
print("\nColumns/Datasets information")
display(All_Metadata_col.head())
```

```py
df_real = cohort1.merge_data("data_matrix")
print("\nData Matrix")
display(df_real.head())
```

## Tutorial Notebooks

1. [Creating Multiple Cohorts in TCGA](https://github.com/ElucidataInc/polly-python/blob/main/Enrich/cohort_creation_demo.ipynb)

2. [Proteomics Data Analysis in TCGA using Cohorts](https://github.com/ElucidataInc/polly-python/blob/main/Analyse/consumption_starter_notebooks/Proteomics_polly_python_TCGA.ipynb)
