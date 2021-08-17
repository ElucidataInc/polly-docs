# Cell-type annotation of single cell transcriptomics data using marker based methods

## Introduction

In the data processing protocols of scRNA-seq experiments, cell type identification is a vital step for subsequent analysis.
Mainly, two types of strategies have been reported, e.g., cell-based and cluster-based annotation.
The major challenge of cell-based strategy lies in the determination of cell types on each cluster as multiple cells with different types are present in one cluster.

Cluster-based methods use the pre-computed clusters to identify markers within each cluster. The markers are then use for cell type identification from reference databases which contain a plethora of manually curated marker to cell-type mappings.


### Import functions


```sos
suppressMessages(library(Seurat))
suppressMessages(library(anndata))
suppressMessages(library(dplyr))
suppressMessages(library(clustermole))
```

Defining custom functions to use


```sos
get_seurat_obj = function(adata) {
    
    # Get the expression matrix
    exprs <- t(as.matrix(adata$X))
    colnames(exprs) <- adata$obs_names
    rownames(exprs) <- adata$var_names
    # Create the Seurat object
    seurat <- CreateSeuratObject(exprs)
    # Set the expression assay
    seurat <- SetAssayData(seurat, "data", exprs)
    # Add observation metadata
    seurat <- AddMetaData(seurat, adata$obs)
    # Add embedding
    embedding_2 <- as.list(adata$obsm["X_tsne"])
    rownames(embedding_2$X_tsne) <- adata$obs_names
    colnames(embedding_2$X_tsne) <- c("tsne_1", "tsne_2")
    seurat[["tsne"]] <- CreateDimReducObject(embedding_2$X_tsne, key = "tsne_")
    
    return(seurat)
    
}
```

# 1. Querying Biomedical Molecular data

### 1.1 Install Polly Python


```sos
!pip3 install polly-python --user
```

    Looking in indexes: https://pypi.org/simple, http://54.245.179.143:80/
    Collecting polly-python
      Downloading https://files.pythonhosted.org/packages/a5/5c/3a3d40cdd4d61e99e7da4a0c2899a6d5e81c093667abf25aba2f157d18db/polly_python-0.0.4-py3-none-any.whl
    Requirement already satisfied: pandas in /usr/local/lib/python3.6/dist-packages (from polly-python) (1.1.5)
    Requirement already satisfied: requests in /usr/local/lib/python3.6/dist-packages (from polly-python) (2.25.1)
    Requirement already satisfied: certifi in /usr/local/lib/python3.6/dist-packages (from polly-python) (2020.6.20)
    Requirement already satisfied: idna in /usr/local/lib/python3.6/dist-packages (from polly-python) (2.8)
    Collecting postpy2 (from polly-python)
      Downloading https://files.pythonhosted.org/packages/7e/25/d98f5ea87d937bf44963be78b504dbb3ddea1c6ddfcc0985b02d062b41dd/postpy2-0.0.6-py3-none-any.whl
    Collecting python-magic (from polly-python)
      Downloading https://files.pythonhosted.org/packages/d3/99/c89223c6547df268596899334ee77b3051f606077317023617b1c43162fb/python_magic-0.4.24-py2.py3-none-any.whl
    Requirement already satisfied: urllib3 in /usr/local/lib/python3.6/dist-packages (from polly-python) (1.26.3)
    Requirement already satisfied: chardet in /usr/local/lib/python3.6/dist-packages (from polly-python) (3.0.4)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas->polly-python) (2019.3)
    Requirement already satisfied: python-dateutil>=2.7.3 in /usr/local/lib/python3.6/dist-packages (from pandas->polly-python) (2.8.0)
    Requirement already satisfied: numpy>=1.15.4 in /usr/local/lib/python3.6/dist-packages (from pandas->polly-python) (1.19.5)
    Requirement already satisfied: six>=1.5 in /usr/local/lib/python3.6/dist-packages (from python-dateutil>=2.7.3->pandas->polly-python) (1.14.0)
    Installing collected packages: python-magic, postpy2, polly-python
    Successfully installed polly-python-0.0.4 postpy2-0.0.6 python-magic-0.4.24
    [33mWARNING: You are using pip version 19.2.3, however version 21.1.3 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.[0m



```sos
from polly.omixatlas import OmixAtlas
import pandas as pd
import os

LiverOmixAltasID = "1615965444377"
```

### 1.2 Connect to Polly OmixAtlas

**Get Authentication Tokens**


```sos
REFRESH_TOKEN="eyJjdHkiOiJKV1QiLCJlbmMiOiJBMjU2R0NNIiwiYWxnIjoiUlNBLU9BRVAifQ.PjS7n0ExVq2ciJnp_uHu5SW2f1Bx1a56aDQNXKgYbVr62DaFlPbcG-ETZ4X6CH2kefjUc4U0h2hzri_lMewpobm0-NKbZJUk0-XBpcNYW0aqEuSdDN_lgnNXIKjTFWHilPQfnweqHKwcmxC2FBOvaj6w3rk0Um3TY2L7H5cEppGe0pqdyc2nGYPNIJ6O9zNMtXxVPkb_TqUVSye48wG7hSNE7nafXZIirZlhvMnuSdYOycKqti9XJbmSB2DB7oMfeqJyiGXywA6c5HSpaxfSO_Mq09XKaXrCnA8ZvRl-v-OEmfhjJCabJeRPSLOApMH7TUhYLc7CME_VMekK0QpWgQ.y4irbtLWSaMtWwdk.Wj7Y2YioY7m9qaoRc934IfZ6CpQ5UXdttWb8onXpsHKqWEc-KTRufdDZPrvOvJc3bJ9InULISS0GAtQzb2wZ_t-1JOKkujm6pl1hZZKl7ynm7cfx4xE23XX8_Qlp6eZsIHn1RQWm9NZq2SzjU5NZWYYLhiSwuO_yElIuPouUanNJqRoXqFK4oKaIAKKHJvPC4iH765qmTiKBJ-StMfvKSyED4kskwA1GJnrdwMd8aP-BmuXvPnU0MzI91wkS4sQIIZcXnbG23INUGJmt_5t3xzCJjgDt-tbRFl0w2949726HJkDPqrVlsOs5uiXbtLNu1UQt94HDscqbIl0BgOAYwAZEvEFa5hRnAGElo5g7AM2YkpX7pVcz1Ne4OegCwvYrDUaH83xoUguNHrC13uQxtjaGqJgNGSDKWuadZdAA1J9YU0uMI8Qngsb-8Ed2Pgr0XArpO1sQhamlDrkjDA6qXXMXT3rM7dEdiKvrtIyWKDim6wcGesVO8FMHXl_dsT9q_991ePVf5c0N82rUohh-0qW-uUGZ350uk0n_KxHdjFDPYr_GZIzmIlxxMegGSgfqe5Z0KbwzWp-a11_1wpO5ts3YQOYmUofzemCna3Psv0ooYNf6qF9KtIGG71rgOJXtxohAwijc1uo7B0MJVAydP3nIFz1GSyYh4dy6WwToeDMKVAvnbAR5MSDVC0fI1orOVuhH0JAZsBq0lV7n5dSC5GmWjAuJkRrHKyOS17LcrgLsmarDVhThXs-84hP2pvyI0TBIHrLfeflstgauT0FEu_z4lzyfGyvkHjJw7V6pyjrRWgUzbP6NNbNVG3njMLPG0Gjm0yaieIzWBtRctVVOZjsuhQTkt7C3nbJPSqC3LApLt_bRDrDIqrLJkrwBQLX9CRsVM72DnSHBQu7nYPP6v14F45ZGtToORWNJ3vf5v0dJqRWWx6YBVXjgx77wylAM98PH4E4LcKysjvmranma-8WQ8l9yIRz7HmnIGbYGIE6xU7KSRvPtgRdVKoD9zzqk5hYn4ncquweDi6o-9sRxs6dM7HMslznCi-mzrdfH_xw43LRTAH-9QWSAAuAnUNS1oShVhHlZhzNjVxHHKq9rD_Ot5bA1CD8Mt-h_P9Hf_S65-3fKCSSrp8xXCe0EkfqYJfYO2VsjiRBPiWFKXEU_4ywNv9TQiszqQaDVdH8TV4PMevWdMPj-tLRNh5D81YRlh7wF8eq8Tj9DFVI5zbMbnQsgHmVzxXJ8G-SsJW41i1P5_BZcwDMaOFWAvdLtUsIHEvJjo8Ut6SRMlx7Pmn0TcBEFzB8oPwHrm9KSz26NfoEnpRidn84sCfs8tDM.piZrXCf1tDWIFBt0RrTPYg"
repo_client = OmixAtlas(REFRESH_TOKEN)
```

**See a list of available Omix Atlases**


```sos
pd.DataFrame.from_dict(repo_client.get_all_omixatlas()['data'])
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>repo_name</th>
      <th>repo_id</th>
      <th>indexes</th>
      <th>dataset_count</th>
      <th>disease_count</th>
      <th>diseases</th>
      <th>organism_count</th>
      <th>organisms</th>
      <th>sources</th>
      <th>datatypes</th>
      <th>sample_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>liveromix_atlas</td>
      <td>1615965444377</td>
      <td>{'gct_metadata': 'liveromix_atlas_gct_metadata...</td>
      <td>6760</td>
      <td>761</td>
      <td>[normal, carcinoma, hepatocellular, obesity, n...</td>
      <td>22</td>
      <td>[homo sapiens, mus musculus, rattus norvegicus...</td>
      <td>[geo, lincs, tcga, metabolomics workbench, met...</td>
      <td>[transcriptomics, mutation, metabolomics, sing...</td>
      <td>1738079</td>
    </tr>
  </tbody>
</table>
</div>



### 1.3 Querying Single Cell data from Liver OmixAtlas

The publication provides the GEO dataset id that was used for the analysis. Using the GSE and GPL id the GEO datasets are queried on the Liver OmixAtlas.  
The data is available in a structured anndata object in `h5ad` format 


```sos
sc_normal_query = """select * FROM liveromix_atlas_files 
        WHERE dataset_id = 'GSE115469_GPL16791'
        AND kw_data_type = 'Single cell'"""

repo_client.query_metadata(sc_normal_query)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>dataset_id</th>
      <th>dataset_source</th>
      <th>disease</th>
      <th>description</th>
      <th>organism</th>
      <th>tissue</th>
      <th>file_type</th>
      <th>authors</th>
      <th>year</th>
      <th>publication</th>
      <th>kw_data_type</th>
      <th>kw_cell_type</th>
      <th>platform</th>
      <th>kw_repo</th>
      <th>kw_package</th>
      <th>kw_key</th>
      <th>kw_bucket</th>
      <th>kw_filetype</th>
      <th>kw_region</th>
      <th>kw_location</th>
      <th>kw_timestamp</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GSE115469_GPL16791</td>
      <td>GEO</td>
      <td>Normal</td>
      <td>Single cell RNA sequencing of human liver reve...</td>
      <td>Homo sapiens</td>
      <td>[liver]</td>
      <td>h5ad</td>
      <td>Sonya A. MacParland, Ian D McGilvray</td>
      <td>2018</td>
      <td>https://www.nature.com/articles/s41467-018-063...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>10x</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE115469_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626622850659</td>
    </tr>
  </tbody>
</table>
</div>




```sos
dataset_id = list(repo_client.query_metadata(sc_normal_query)['dataset_id'])
dataset_id
```




    ['GSE115469_GPL16791']



**Download data from liverOmix atlas**


```sos
file_link = repo_client.download_data(LiverOmixAltasID, dataset_id[0])['data']
file_link
```




    'https://discover-prod-datalake-v1.s3.amazonaws.com/liver_atlas/data/SingleCell/GSE115469_GPL16791.h5ad?AWSAccessKeyId=ASIAVRYB5UBIA5FW2NFB&Signature=wQeC6muM1IPeU2RQGOWf7z2gaMU%3D&x-amz-security-token=IQoJb3JpZ2luX2VjENf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQDQlEn4h5AZ3rKsV3RrcKVNwYuMCpB8CW4v0rYUyyCUUwIhALNparKajiikPLlUbV7tVAXgdOqRsCMy%2BEpBKeg8QtflKqoCCND%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQAhoMMzgxNzE5MjU3MTY4Igz90u%2BAalZnFYSf59wq%2FgErXZhA1mS24Pzmxiujn5CDBTadtH%2BpJ%2BC4gLUAjzk3eAFnQuKfZfMN6MNdTN81pvl3g3p9Cr1LLeiPF%2FxHNEdltiie79cRDd3dZF0JMUqbWUY7y6%2Fy9It2YluhF2JGAGf8PqPT9zadygDz%2FYRLzRryb6Xcka0t%2FD6ST6P%2BZ6UlOxplcslONUukwhWgaIdPhl4cwUSP7LCL41fJAxUCiEb8%2FugQMMPvLu62qFGfdQpG7bUFddfXqSuE3mST0EAlY5b64c46YjgbadhJd52e9uIowubqpMxVsCOt0VCpLFE54ITnrMnv9t0Qp4KsNJHBu%2F4%2Bq9%2BvfxX%2BFkitnHj7oTCvtOSHBjqZAVP5nfMsDhkbIXsOcKWxas6yjg5q9eZDOecxXzT%2BR0iRAcMLh6fVhYciG327iZOaq0NQ87LDxkgT3zY9f4g0EN4DSj3jE5O28lrSPyOzbbuB9rhzJuu1A3uaap00bTErXKquKUFM4cPgCpbADzIoQQLuzoX6JSOXWySZ%2BN%2FX0MMFa4eWiZAUZyDkoet9oY4%2B%2BcdkndKPQIrceg%3D%3D&Expires=1626941537'




```sos
wget -O GSE115469_GPL16791.h5ad 'https://discover-prod-datalake-v1.s3.amazonaws.com/liver_atlas/data/SingleCell/GSE115469_GPL16791.h5ad?AWSAccessKeyId=ASIAVRYB5UBIA5FW2NFB&Signature=wQeC6muM1IPeU2RQGOWf7z2gaMU%3D&x-amz-security-token=IQoJb3JpZ2luX2VjENf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQDQlEn4h5AZ3rKsV3RrcKVNwYuMCpB8CW4v0rYUyyCUUwIhALNparKajiikPLlUbV7tVAXgdOqRsCMy%2BEpBKeg8QtflKqoCCND%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQAhoMMzgxNzE5MjU3MTY4Igz90u%2BAalZnFYSf59wq%2FgErXZhA1mS24Pzmxiujn5CDBTadtH%2BpJ%2BC4gLUAjzk3eAFnQuKfZfMN6MNdTN81pvl3g3p9Cr1LLeiPF%2FxHNEdltiie79cRDd3dZF0JMUqbWUY7y6%2Fy9It2YluhF2JGAGf8PqPT9zadygDz%2FYRLzRryb6Xcka0t%2FD6ST6P%2BZ6UlOxplcslONUukwhWgaIdPhl4cwUSP7LCL41fJAxUCiEb8%2FugQMMPvLu62qFGfdQpG7bUFddfXqSuE3mST0EAlY5b64c46YjgbadhJd52e9uIowubqpMxVsCOt0VCpLFE54ITnrMnv9t0Qp4KsNJHBu%2F4%2Bq9%2BvfxX%2BFkitnHj7oTCvtOSHBjqZAVP5nfMsDhkbIXsOcKWxas6yjg5q9eZDOecxXzT%2BR0iRAcMLh6fVhYciG327iZOaq0NQ87LDxkgT3zY9f4g0EN4DSj3jE5O28lrSPyOzbbuB9rhzJuu1A3uaap00bTErXKquKUFM4cPgCpbADzIoQQLuzoX6JSOXWySZ%2BN%2FX0MMFa4eWiZAUZyDkoet9oY4%2B%2BcdkndKPQIrceg%3D%3D&Expires=1626941537'
```

    --2021-07-22 07:12:54--  https://discover-prod-datalake-v1.s3.amazonaws.com/liver_atlas/data/SingleCell/GSE115469_GPL16791.h5ad?AWSAccessKeyId=ASIAVRYB5UBIA5FW2NFB&Signature=wQeC6muM1IPeU2RQGOWf7z2gaMU%3D&x-amz-security-token=IQoJb3JpZ2luX2VjENf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQDQlEn4h5AZ3rKsV3RrcKVNwYuMCpB8CW4v0rYUyyCUUwIhALNparKajiikPLlUbV7tVAXgdOqRsCMy%2BEpBKeg8QtflKqoCCND%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQAhoMMzgxNzE5MjU3MTY4Igz90u%2BAalZnFYSf59wq%2FgErXZhA1mS24Pzmxiujn5CDBTadtH%2BpJ%2BC4gLUAjzk3eAFnQuKfZfMN6MNdTN81pvl3g3p9Cr1LLeiPF%2FxHNEdltiie79cRDd3dZF0JMUqbWUY7y6%2Fy9It2YluhF2JGAGf8PqPT9zadygDz%2FYRLzRryb6Xcka0t%2FD6ST6P%2BZ6UlOxplcslONUukwhWgaIdPhl4cwUSP7LCL41fJAxUCiEb8%2FugQMMPvLu62qFGfdQpG7bUFddfXqSuE3mST0EAlY5b64c46YjgbadhJd52e9uIowubqpMxVsCOt0VCpLFE54ITnrMnv9t0Qp4KsNJHBu%2F4%2Bq9%2BvfxX%2BFkitnHj7oTCvtOSHBjqZAVP5nfMsDhkbIXsOcKWxas6yjg5q9eZDOecxXzT%2BR0iRAcMLh6fVhYciG327iZOaq0NQ87LDxkgT3zY9f4g0EN4DSj3jE5O28lrSPyOzbbuB9rhzJuu1A3uaap00bTErXKquKUFM4cPgCpbADzIoQQLuzoX6JSOXWySZ%2BN%2FX0MMFa4eWiZAUZyDkoet9oY4%2B%2BcdkndKPQIrceg%3D%3D&Expires=1626941537
    Resolving discover-prod-datalake-v1.s3.amazonaws.com (discover-prod-datalake-v1.s3.amazonaws.com)... 52.218.216.170
    Connecting to discover-prod-datalake-v1.s3.amazonaws.com (discover-prod-datalake-v1.s3.amazonaws.com)|52.218.216.170|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 186223226 (178M) [binary/octet-stream]
    Saving to: ‘GSE115469_GPL16791.h5ad’
    
    GSE115469_GPL16791. 100%[===================>] 177.60M  30.7MB/s    in 6.3s    
    
    2021-07-22 07:13:01 (28.0 MB/s) - ‘GSE115469_GPL16791.h5ad’ saved [186223226/186223226]
    


## 2. Data

The single-cell dataset used is a map of the cellular landscape of the human liver using single-cell RNA sequencing. The dataset consists of transcriptional profiles of 8444 parenchymal and non-parenchymal cells obtained from the fractionation of fresh hepatic tissue from five human livers.  

The original publication shows 12 distinct cell type populations identified using gene expression patterns, flow cytometry, and immunohistochemical examinations. However, the dataset deposited on Gene Expression Omnibus (GEO), did not contain the experimentally identified cell type annotations.


```sos
options(warn=-1)

adata.pred = read_h5ad('GSE115469_GPL16791.h5ad')

seurat = get_seurat_obj(adata.pred)

DimPlot(object = seurat, reduction = "tsne", group.by = "clusters", cols = DiscretePalette(20), label = TRUE)
```


![png](output_24_0.png)


## 3. Finding cluster biomarkers

A one vs all comparison is performed between clusters to identify marker genes.


```sos
# find markers for every cluster compared to all remaining cells, report only the positive ones
Idents(object = seurat) <- seurat@meta.data$clusters
liver.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```


```sos
head(liver.markers)
```


<table class="dataframe">
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>HSD11B1</th><td>0</td><td>0.7473840</td><td>0.896</td><td>0.207</td><td>0</td><td>1</td><td>HSD11B1</td></tr>
	<tr><th scope=row>APOM</th><td>0</td><td>0.6909939</td><td>0.913</td><td>0.204</td><td>0</td><td>1</td><td>APOM   </td></tr>
	<tr><th scope=row>PON3</th><td>0</td><td>0.5102462</td><td>0.858</td><td>0.177</td><td>0</td><td>1</td><td>PON3   </td></tr>
	<tr><th scope=row>TTC36</th><td>0</td><td>0.4761022</td><td>0.827</td><td>0.193</td><td>0</td><td>1</td><td>TTC36  </td></tr>
	<tr><th scope=row>F10</th><td>0</td><td>0.4080488</td><td>0.583</td><td>0.104</td><td>0</td><td>1</td><td>F10    </td></tr>
	<tr><th scope=row>BCHE</th><td>0</td><td>0.3824872</td><td>0.749</td><td>0.156</td><td>0</td><td>1</td><td>BCHE   </td></tr>
</tbody>
</table>



**Plotting a heatmap showing the top 5 clusters containing largest number of cells**


```sos
#Scaling data
seurat <- ScaleData(seurat)

Idents(object = seurat) <- seurat@meta.data$clusters

liver.markers.top <- liver.markers %>% group_by(cluster) %>%  top_n(n = 2 ,wt = avg_log2FC)
liver.markers.subset <- liver.markers %>% group_by(cluster) %>% filter(cluster == c(1,2,3,4,5)) %>% top_n(n = 2 ,wt = avg_log2FC)

seurat.subset <- subset(seurat, idents = c('1','2','3','4','5'))
DoHeatmap(seurat.subset, features = liver.markers.subset$gene)
```

    Centering and scaling data matrix
    



![png](output_30_1.png)


## 4 . Automated cell type curation using biomarker genes

Marker gene database-based annotation takes advantage of cell type atlases. Literature- and scRNA-seq analysis-derived markers have been assembled into reference cell type hierarchies and marker lists. In this approach, basic scoring systems are used to ascribe cell types at the cluster level in the query dataset.


```sos
Type <- c()
clusters <- c()
for (i in seq(1,20, by=1)) {
    liver.markers.filt <- as.data.frame(liver.markers %>% filter(cluster == as.character(i)) %>% top_n(n = 25, wt = avg_log2FC))
    my_overlaps <- clustermole_overlaps(genes = liver.markers.filt$gene, species = "hs")
    type <- as.character(my_overlaps[1,5])
    Type <- c(Type, type)
    clusters <- c(clusters, i)
}


clust_map <- as.data.frame(cbind(clusters,Type))
head(clust_map)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>clusters</th><th scope=col>Type</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>1</td><td>Hepatocytes                                        </td></tr>
	<tr><th scope=row>2</th><td>2</td><td>DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_CD8_T_CELLS</td></tr>
	<tr><th scope=row>3</th><td>3</td><td>Hepatocytes                                        </td></tr>
	<tr><th scope=row>4</th><td>4</td><td>AIZARANI_LIVER_C25_KUPFFER_CELLS_4                 </td></tr>
	<tr><th scope=row>5</th><td>5</td><td>Hepatocytes                                        </td></tr>
	<tr><th scope=row>6</th><td>6</td><td>AIZARANI_LIVER_C11_HEPATOCYTES_1                   </td></tr>
</tbody>
</table>



### Curating cell types from cluster based identifications

The cell type annotations predicted lacked harmonisation in cell names, since the algorithm used a marker-gene database search which queries markers from multiple databases each with database specific vocabulary. Therefore, in order for the cell types to cluster, markers corresponding to same cell type but different annotation must be curated to follow a common vocabulary


```sos
#The cell type names must be cleaned in order for them to be clustered together
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('NK_NKT_CELLS','NK Cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('B cell\\D+','B cell',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_CD71\\D+','Erythroid cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_NK_CELLS','NK Cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('KUPFFER_CELLS','Kupffer Cells',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('MVECS','Endothelial cells',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_CD8_T_CELLS','CD8 T cell',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_CD4_T_CELLS','CD4 T cell',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('AIZARANI_LIVER_C\\d+_','',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('_\\d+$','',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('HEPATOCYTES','Hepatocytes',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\bPlasma cell\\b','Plasma cells',x)) 
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_PLASMA_CELLS','Plasma cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_PLASMA_CELL','Plasma cells',x))                               
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('EPCAM_POS_\\D+','EPCAM+',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_B_CELL','B cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('CUI_DEVELOPING_HEART_C3_FIBROBLAST_LIKE_CELL','Others',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('Erythrocytes_\\D+','Erythrocytes',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('FAN_EMBRYONIC_CTX_NSC','Others',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('Lymph_endothel_cell','Other endothelial cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('\\D+_Epithelial_cells','Epithelial cells',x))
clust_map$Type <- sapply(clust_map$Type, function(x) gsub('HAY_BONE_MARROW_PRO_B','Others',x))

                               
cell.type_clust <- aggregate(clusters ~ Type, clust_map, paste, collapse=",")
cell.type_clust
```


<table class="dataframe">
<caption>A data.frame: 12 × 2</caption>
<thead>
	<tr><th scope=col>Type</th><th scope=col>clusters</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>B cell                </td><td>16           </td></tr>
	<tr><td>CD8 T cell            </td><td>2            </td></tr>
	<tr><td>Cholangiocytes        </td><td>17           </td></tr>
	<tr><td>Endothelial cells     </td><td>11,13        </td></tr>
	<tr><td>Erythroid cells       </td><td>19           </td></tr>
	<tr><td>Gamma delta T cells   </td><td>18           </td></tr>
	<tr><td>Hepatic stellate cells</td><td>20           </td></tr>
	<tr><td>Hepatocytes           </td><td>1,3,5,6,14,15</td></tr>
	<tr><td>Kupffer Cells         </td><td>4,10         </td></tr>
	<tr><td>LSECS                 </td><td>12           </td></tr>
	<tr><td>NK Cells              </td><td>8,9          </td></tr>
	<tr><td>Plasma cells          </td><td>7            </td></tr>
</tbody>
</table>




```sos
rows <- rownames(seurat@meta.data)
seurat@meta.data = merge(seurat@meta.data, clust_map, by = 'clusters', all.x = TRUE)
rownames(seurat@meta.data) <- rows
```


```sos
DimPlot(object = seurat, reduction = "tsne", group.by = "Type", cols = DiscretePalette(12))
```


![png](output_36_0.png)


## Appendix

### Package installation


```sos
!sudo R -e 'BiocManager::install(c("GSEABase", "GSVA", "singscore"))'
```


```sos
!sudo R -e 'install.packages(c("anndata", "clustermole"), repos = "https://cloud.r-project.org/")'
```
