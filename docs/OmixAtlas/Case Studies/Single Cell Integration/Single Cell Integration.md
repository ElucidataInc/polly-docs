# Integrative single-cell transcriptome analysis reveals a subpopulation of fibroblasts associated with favorable prognosis of liver cancer patients

#### Authors *H. Wang, C. Feng, M. Lu*

## Introduction

An integrative single cell transcriptomic data was performed between normal and tumor patients from the liver tissue. Here we show the power of integration-ready biomedical molecular data to identify markers associated with fibroblast subpopulation. 

Normal and tumor single cell datasets are shown to *readily integrate* due to the standardised data structure on Liver OmixAtlas. 

In order to evaluate the clinical significance of those markers, we utilise curated transcriptomics data for liver hepatocellular carcinoma to perform survival analysis in cancer patients. The overall analysis showed that overexpression of **SPARCL1** is associated with favourable prognosis of liver cancer patients. 

### Import functions

#### Please follow the instructions below

**1. Check the appendix section below to install necessary packages**  
**2. Refresh the kernel after installation is complete**  


```sos
import re
from cmapPy.pandasGEXpress.parse_gct import parse
from matplotlib.pyplot import rc_context
import itertools
import pandas as pd
import numpy as np
import GEOparse
from sh import gunzip
import wget
import urllib.request
import os
```

Defining custom functions to use


```sos
def get_data_geo(accession):

    if not os.path.isfile('./metadata/'+accession+'_family.soft.gz'): 
    
        GEOparse.get_GEO(geo=accession, destdir="./metadata/", silent=True)  
        gunzip('./metadata/'+accession+'_family.soft.gz')

    os.makedirs(accession, exist_ok = True)

    with open('./metadata/'+accession+'_family.soft', 'rt') as f:
        for line in f:
            if('!Series_supplementary_file' in line):
                file = str(line.partition("=")[2]).strip()
                urllib.request.urlretrieve(file, accession+'/'+file.split('/')[-1])
            else:
                continue
    f.close()
```

# 1. Querying Biomedical Molecular data

## 1.1 On Polly Liver OmixAtlas

### 1.1.1 Install Polly Python


```sos
!pip3 install polly-python --user
```

    Looking in indexes: https://pypi.org/simple, http://54.245.179.143:80/
    Collecting polly-python
      Downloading https://files.pythonhosted.org/packages/a5/5c/3a3d40cdd4d61e99e7da4a0c2899a6d5e81c093667abf25aba2f157d18db/polly_python-0.0.4-py3-none-any.whl
    Requirement already satisfied: idna in /usr/local/lib/python3.6/dist-packages (from polly-python) (2.8)
    Collecting python-magic (from polly-python)
      Downloading https://files.pythonhosted.org/packages/d3/99/c89223c6547df268596899334ee77b3051f606077317023617b1c43162fb/python_magic-0.4.24-py2.py3-none-any.whl
    Collecting postpy2 (from polly-python)
      Downloading https://files.pythonhosted.org/packages/7e/25/d98f5ea87d937bf44963be78b504dbb3ddea1c6ddfcc0985b02d062b41dd/postpy2-0.0.6-py3-none-any.whl
    Requirement already satisfied: pandas in /usr/local/lib/python3.6/dist-packages (from polly-python) (1.1.5)
    Requirement already satisfied: certifi in /usr/local/lib/python3.6/dist-packages (from polly-python) (2020.6.20)
    Requirement already satisfied: chardet in /usr/local/lib/python3.6/dist-packages (from polly-python) (3.0.4)
    Requirement already satisfied: urllib3 in /usr/local/lib/python3.6/dist-packages (from polly-python) (1.26.3)
    Requirement already satisfied: requests in /usr/local/lib/python3.6/dist-packages (from polly-python) (2.25.1)
    Requirement already satisfied: python-dateutil>=2.7.3 in /usr/local/lib/python3.6/dist-packages (from pandas->polly-python) (2.8.0)
    Requirement already satisfied: numpy>=1.15.4 in /usr/local/lib/python3.6/dist-packages (from pandas->polly-python) (1.19.5)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas->polly-python) (2019.3)
    Requirement already satisfied: six>=1.5 in /usr/local/lib/python3.6/dist-packages (from python-dateutil>=2.7.3->pandas->polly-python) (1.14.0)
    Installing collected packages: python-magic, postpy2, polly-python
    Successfully installed polly-python-0.0.4 postpy2-0.0.6 python-magic-0.4.24
    [33mWARNING: You are using pip version 19.2.3, however version 21.1.3 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.[0m



```sos
from polly.omixatlas import OmixAtlas
```

### 1.1.2 Connect to Polly OmixAtlas  

**Get Authentication Tokens**


```sos
#REFRESH_TOKEN = "ENTER REFRESH TOKEN WITHIN THE DOUBLE QUOTES"
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



### 1.1.3 Get Liver OmixAtlas Summary


```sos
LiverOmixAltasID = "1615965444377"
LiverOmixAtlas_summary = repo_client.omixatlas_summary(LiverOmixAltasID)
pd.DataFrame.from_dict(LiverOmixAtlas_summary, orient="index")
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
      <th>data</th>
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




```sos
#Data sources in Liver OmixAtlas
LiverOmixAtlas_summary['data'].get('sources')
```




    ['geo',
     'lincs',
     'tcga',
     'metabolomics workbench',
     'metabolights',
     'ccle',
     'gene expression omnibus (geo)',
     'cptac',
     'depmap',
     'human protein atlas']




```sos
#Data Types in Liver OmixAtlas
LiverOmixAtlas_summary['data'].get('datatypes')
```




    ['transcriptomics',
     'mutation',
     'metabolomics',
     'single cell',
     'proteomics',
     'lipidomics',
     'mirna',
     'drug screens',
     'gene dependency',
     'gene effect']




```sos
LiverOmixAtlas_summary['data'].get('indexes')
```




    {'gct_metadata': 'liveromix_atlas_gct_metadata',
     'h5ad_metadata': 'liveromix_atlas_h5ad_metadata',
     'csv': 'liveromix_atlas_csv',
     'files': 'liveromix_atlas_files',
     'json': 'liveromix_atlas_json',
     'ipynb': 'liveromix_atlas_ipynb',
     'gct_data': 'liveromix_atlas_gct_data',
     'h5ad_data': 'liveromix_atlas_h5ad_data'}



### 1.1.4 Querying Single Cell RNASeq data on Liver OmixAtlas

**All data in Liver Atlas are structured with metadata fields tagged with ontologies and controlled vocabularies, which simplifies finding relevant datasets**

To access, filter, and search the metadata in Liver OmixAtlas, the function mentioned below can be used: `query_metadata` (“query written in SQL”)

For integrative single cell analysis, the authors used single cell transcriptomics data from:  
1. 5 normal patients [GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469)  
2. 12 cancer patients [GSE125449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449)  

The study investigated a total of 8439 cells from normal patients and and 9946 cells from liver cancer tissues.
  
For querying the Single Cell RNASeq datasets we utilised the GEO Accession ID that was made available in the publication.   

**Fetching liver tumor dataset**


```sos
sc_tumor_query = """select * FROM liveromix_atlas_files 
        WHERE kw_data_type = 'Single cell'
        AND organism = 'Homo sapiens'
        AND (disease = 'Carcinoma, Hepatocellular'
            OR disease = 'Tumor')
        """

sc_tumor_query_df = repo_client.query_metadata(sc_tumor_query)
```

**In this case the query returned the datasets that are required, but requires additional filtering**


```sos
sc_tumor_query_df
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
      <th>file_type</th>
      <th>tissue</th>
      <th>disease</th>
      <th>dataset_id</th>
      <th>organism</th>
      <th>dataset_source</th>
      <th>platform</th>
      <th>description</th>
      <th>kw_data_type</th>
      <th>kw_cell_type</th>
      <th>curation_version</th>
      <th>source_process</th>
      <th>publication</th>
      <th>kw_cell_line</th>
      <th>kw_drug</th>
      <th>kw_gene</th>
      <th>kw_repo</th>
      <th>kw_package</th>
      <th>kw_key</th>
      <th>kw_bucket</th>
      <th>kw_filetype</th>
      <th>kw_region</th>
      <th>kw_location</th>
      <th>kw_timestamp</th>
      <th>summary</th>
      <th>overall_design</th>
      <th>pubmed_id</th>
      <th>total_num_samples</th>
      <th>total_num_cells</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>h5ad</td>
      <td>[liver, peripheral blood]</td>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>GSE98638_GPL16791</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL16791</td>
      <td>[Landscape of infiltrating T cells in liver ca...</td>
      <td>Single cell</td>
      <td>[helper T cell, regulatory T cell]</td>
      <td>g2</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[CD8A, CD4, IL2RA, TTR, PTCHD3, PTH]</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE98638_GPL16791....</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626796397637</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>[Cholangiocarcinoma, Carcinoma, Hepatocellular]</td>
      <td>GSE125449_GPL20301</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>GPL20301</td>
      <td>Tumor cell biodiversity drives microenvironmen...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE125449_GPL20301...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795852462</td>
      <td>Single-cell transcriptome profiling of liver c...</td>
      <td>[A total of 19 tumors were profiled. Set 1 con...</td>
      <td>31588021</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>Tumor</td>
      <td>GSE125449_GPL18573</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>GPL18573</td>
      <td>Tumor cell biodiversity drives microenvironmen...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE125449_GPL18573...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626622880810</td>
      <td>Single-cell transcriptome profiling of liver c...</td>
      <td>[A total of 19 tumors were profiled. Set 1 con...</td>
      <td>31588021</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>h5ad</td>
      <td>[blood]</td>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>GSE107747_GPL16791</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL16791</td>
      <td>Single-cell RNA sequencing of circulating tumo...</td>
      <td>Single cell</td>
      <td>[neoplastic cell]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE107747_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795662655</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2</td>
      <td>10184</td>
    </tr>
    <tr>
      <th>4</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>GSE112271_GPL16791</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL16791</td>
      <td>Single-cell RNA of multiregional sampling in h...</td>
      <td>Single cell</td>
      <td>[None]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE112271_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795689768</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>7</td>
      <td>31279</td>
    </tr>
  </tbody>
</table>
</div>



**We perform a partial match using the GSE ID provided in the publication**


```sos
sc_tumor_query_df_filt = sc_tumor_query_df[sc_tumor_query_df['dataset_id'].str.contains('GSE125449')]
sc_tumor_query_df_filt
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
      <th>file_type</th>
      <th>tissue</th>
      <th>disease</th>
      <th>dataset_id</th>
      <th>organism</th>
      <th>dataset_source</th>
      <th>platform</th>
      <th>description</th>
      <th>kw_data_type</th>
      <th>kw_cell_type</th>
      <th>curation_version</th>
      <th>source_process</th>
      <th>publication</th>
      <th>kw_cell_line</th>
      <th>kw_drug</th>
      <th>kw_gene</th>
      <th>kw_repo</th>
      <th>kw_package</th>
      <th>kw_key</th>
      <th>kw_bucket</th>
      <th>kw_filetype</th>
      <th>kw_region</th>
      <th>kw_location</th>
      <th>kw_timestamp</th>
      <th>summary</th>
      <th>overall_design</th>
      <th>pubmed_id</th>
      <th>total_num_samples</th>
      <th>total_num_cells</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>[Cholangiocarcinoma, Carcinoma, Hepatocellular]</td>
      <td>GSE125449_GPL20301</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>GPL20301</td>
      <td>Tumor cell biodiversity drives microenvironmen...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE125449_GPL20301...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795852462</td>
      <td>Single-cell transcriptome profiling of liver c...</td>
      <td>[A total of 19 tumors were profiled. Set 1 con...</td>
      <td>31588021</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>Tumor</td>
      <td>GSE125449_GPL18573</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>GPL18573</td>
      <td>Tumor cell biodiversity drives microenvironmen...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE125449_GPL18573...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626622880810</td>
      <td>Single-cell transcriptome profiling of liver c...</td>
      <td>[A total of 19 tumors were profiled. Set 1 con...</td>
      <td>31588021</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```sos
sc_tumor_id = list(sc_tumor_query_df_filt['dataset_id'])
sc_tumor_id
```




    ['GSE125449_GPL20301', 'GSE125449_GPL18573']



**Download single cell tumor data from OmixAtlas**


```sos
url = repo_client.download_data(LiverOmixAltasID, sc_tumor_id[0]).get('data')
file_name = sc_tumor_id[0]+".h5ad"
os.system(f"wget -O '{file_name}' '{url}'")
```




    0




```sos
url = repo_client.download_data(LiverOmixAltasID, sc_tumor_id[1]).get('data')
file_name = sc_tumor_id[1]+".h5ad"
os.system(f"wget -O '{file_name}' '{url}'")
```




    0



**Fetching normal liver dataset**


```sos
sc_normal_query = """select * FROM liveromix_atlas_files 
        WHERE disease = 'Normal' 
        AND kw_data_type = 'Single cell'
        AND organism = 'Homo sapiens'"""

sc_normal_query_df = repo_client.query_metadata(sc_normal_query)
sc_normal_query_df
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
      <th>file_type</th>
      <th>tissue</th>
      <th>disease</th>
      <th>dataset_id</th>
      <th>organism</th>
      <th>dataset_source</th>
      <th>platform</th>
      <th>description</th>
      <th>kw_data_type</th>
      <th>kw_cell_type</th>
      <th>curation_version</th>
      <th>source_process</th>
      <th>publication</th>
      <th>kw_cell_line</th>
      <th>kw_drug</th>
      <th>kw_gene</th>
      <th>total_num_samples</th>
      <th>total_num_cells</th>
      <th>kw_repo</th>
      <th>kw_package</th>
      <th>kw_key</th>
      <th>kw_bucket</th>
      <th>kw_filetype</th>
      <th>kw_region</th>
      <th>kw_location</th>
      <th>kw_timestamp</th>
      <th>authors</th>
      <th>year</th>
      <th>geo_series_id</th>
      <th>geo_platform_id</th>
      <th>pubmed_ids</th>
      <th>title</th>
      <th>summary</th>
      <th>submission_date</th>
      <th>species</th>
      <th>cell_types</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>h5ad</td>
      <td>[skin, liver, lung]</td>
      <td>[Normal]</td>
      <td>GSE133345_GPL20795</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL20795</td>
      <td>Deciphering human macrophage development at si...</td>
      <td>Single cell</td>
      <td>[tissue-resident macrophage, embryonic cell]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[PTPRC, GYPA]</td>
      <td>9</td>
      <td>1268</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE133345_GPL20795...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626363659950</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>Normal</td>
      <td>GSE115469_GPL16791</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>10x</td>
      <td>Single cell RNA sequencing of human liver reve...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>https://www.nature.com/articles/s41467-018-063...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE115469_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626622850659</td>
      <td>Sonya A. MacParland, Ian D McGilvray</td>
      <td>2018</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>[AIDS Dementia Complex, Normal]</td>
      <td>GSE148796_GPL18573</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL18573</td>
      <td>Identification of pathogenic TRAIL-expressing ...</td>
      <td>Single cell</td>
      <td>[macrophage, T cell, natural killer cell, hema...</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[CD4, RAG2, TNFSF10]</td>
      <td>4</td>
      <td>12129</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE148796_GPL18573...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626796167819</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>h5ad</td>
      <td>[liver, bile]</td>
      <td>[Normal]</td>
      <td>GSE141183_GPL16791</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL16791</td>
      <td>RNA-sequencing with human liver organoids deri...</td>
      <td>Single cell</td>
      <td>[pluripotent stem cell, hepatocyte]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[Bile Acids and Salts, bosentan]</td>
      <td>[None]</td>
      <td>1</td>
      <td>5120</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE141183_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795963571</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>[Normal]</td>
      <td>GSE130073_GPL16791</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL16791</td>
      <td>Single cell RNA-sequencing of iPS cells derive...</td>
      <td>Single cell</td>
      <td>[stem cell, hepatocyte, Kupffer cell, hepatic ...</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[retinoic acid]</td>
      <td>[None]</td>
      <td>1</td>
      <td>3985</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE130073_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795871475</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>5</th>
      <td>h5ad</td>
      <td>[liver, muscle]</td>
      <td>[Normal]</td>
      <td>GSE107552_GPL11154</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL11154</td>
      <td>Mapping human pluripotent stem cell differenti...</td>
      <td>Single cell</td>
      <td>[progenitor cell, stromal cell, pluripotent st...</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[sodium(1+)]</td>
      <td>[None]</td>
      <td>10</td>
      <td>1934</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE107552_GPL11154...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626363174715</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>6</th>
      <td>h5ad</td>
      <td>[skin, spleen, heart, blood, stomach, esophagu...</td>
      <td>[Normal]</td>
      <td>GSE159929_GPL20795</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL20795</td>
      <td>Single-cell transcriptome profiling of an adul...</td>
      <td>Single cell</td>
      <td>[B cell, stromal cell, T cell, fibroblast]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[TRBV20OR9-2, BCR]</td>
      <td>15</td>
      <td>82889</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE159929_GPL20795...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626364686464</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>7</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>Normal</td>
      <td>GSE124395_GPL16791</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>mCelSeq2</td>
      <td>A human liver cell atlas reveals heterogeneity...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>https://www.nature.com/articles/s41586-019-1373-2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE124395_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626622864853</td>
      <td>Nadim Aizarani, Dominic Grün</td>
      <td>2019</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>h5ad</td>
      <td>[spleen, liver, bone marrow, peripheral blood]</td>
      <td>[Normal]</td>
      <td>GSE137539_GPL24676</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL24676</td>
      <td>Single-cell transcriptome profiling reveals ne...</td>
      <td>Single cell</td>
      <td>[granulocyte monocyte progenitor cell, mature ...</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[KIT, CD33]</td>
      <td>3</td>
      <td>29622</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE137539_GPL24676...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626795936357</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>9</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>[Normal]</td>
      <td>GSE150226_GPL19415</td>
      <td>[Homo sapiens, Mus musculus]</td>
      <td>GEO</td>
      <td>GPL19415</td>
      <td>Heterogeneity and chimerism of endothelial cel...</td>
      <td>Single cell</td>
      <td>[endothelial cell]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>6</td>
      <td>3723</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE150226_GPL19415...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626796207901</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>10</th>
      <td>h5ad</td>
      <td>[thymus, liver]</td>
      <td>[Normal]</td>
      <td>GSE133341_GPL20795</td>
      <td>[Homo sapiens]</td>
      <td>GEO</td>
      <td>GPL20795</td>
      <td>Single-cell RNA Sequencing Resolves Spatiotemp...</td>
      <td>Single cell</td>
      <td>[stromal cell, pluripotent stem cell, T cell, ...</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>11</td>
      <td>29722</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE133341_GPL20795...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626363621023</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>11</th>
      <td>h5ad</td>
      <td>[placenta, spleen, thymus, liver, stomach, bon...</td>
      <td>[Normal]</td>
      <td>GSE108097_GPL22245</td>
      <td>[Homo sapiens, Mus musculus]</td>
      <td>GEO</td>
      <td>GPL22245</td>
      <td>Mapping Mouse Cell Atlas by Microwell-seq</td>
      <td>Single cell</td>
      <td>[mesenchymal stem cell]</td>
      <td>g3</td>
      <td>connector</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[ES-E14, CJ7]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>1</td>
      <td>390</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE108097_GPL22245...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626363202288</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>12</th>
      <td>NaN</td>
      <td>[liver, heart]</td>
      <td>[Normal]</td>
      <td>GSE90749_GPL11154</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>NaN</td>
      <td>Induced Pluripotent Stem Cell Differentiation ...</td>
      <td>Single cell</td>
      <td>[pluripotent stem cell, mononuclear cell, hepa...</td>
      <td>g3</td>
      <td>NaN</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[SORT1, PSRC1, CELSR2]</td>
      <td>204</td>
      <td>204</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE90749_GPL11154....</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626796394998</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GSE90749</td>
      <td>GPL11154</td>
      <td>[]</td>
      <td>Induced Pluripotent Stem Cell Differentiation ...</td>
      <td>Genome-wide association studies (GWAS) have hi...</td>
      <td>Dec 01 2016</td>
      <td>Human</td>
      <td>[NA, MSC, Tissue_stem_cells, Fibroblasts, Smoo...</td>
    </tr>
    <tr>
      <th>13</th>
      <td>NaN</td>
      <td>[liver]</td>
      <td>[Normal]</td>
      <td>GSE81252_GPL16791</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>NaN</td>
      <td>Multilineage communication regulates human liv...</td>
      <td>Single cell</td>
      <td>[mesenchymal cell, endothelial cell, mesenchym...</td>
      <td>g3</td>
      <td>NaN</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[VEGFA]</td>
      <td>761</td>
      <td>761</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE81252_GPL16791....</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626796391895</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GSE81252</td>
      <td>GPL16791</td>
      <td>[27669147, 28614297]</td>
      <td>Multilineage communication regulates human liv...</td>
      <td>Conventional 2-D differentiation from pluripot...</td>
      <td>May 09 2016</td>
      <td>Human</td>
      <td>[Hepatocytes, iPS_cells, Endothelial_cells, Sm...</td>
    </tr>
  </tbody>
</table>
</div>



**We perform a partial match using the GSE ID provided in the publication**


```sos
sc_normal_query_df_filt = sc_normal_query_df[sc_normal_query_df['dataset_id'].str.contains('GSE115469')]
sc_normal_query_df_filt
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
      <th>file_type</th>
      <th>tissue</th>
      <th>disease</th>
      <th>dataset_id</th>
      <th>organism</th>
      <th>dataset_source</th>
      <th>platform</th>
      <th>description</th>
      <th>kw_data_type</th>
      <th>kw_cell_type</th>
      <th>curation_version</th>
      <th>source_process</th>
      <th>publication</th>
      <th>kw_cell_line</th>
      <th>kw_drug</th>
      <th>kw_gene</th>
      <th>total_num_samples</th>
      <th>total_num_cells</th>
      <th>kw_repo</th>
      <th>kw_package</th>
      <th>kw_key</th>
      <th>kw_bucket</th>
      <th>kw_filetype</th>
      <th>kw_region</th>
      <th>kw_location</th>
      <th>kw_timestamp</th>
      <th>authors</th>
      <th>year</th>
      <th>geo_series_id</th>
      <th>geo_platform_id</th>
      <th>pubmed_ids</th>
      <th>title</th>
      <th>summary</th>
      <th>submission_date</th>
      <th>species</th>
      <th>cell_types</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>h5ad</td>
      <td>[liver]</td>
      <td>Normal</td>
      <td>GSE115469_GPL16791</td>
      <td>Homo sapiens</td>
      <td>GEO</td>
      <td>10x</td>
      <td>Single cell RNA sequencing of human liver reve...</td>
      <td>Single cell</td>
      <td>[]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>https://www.nature.com/articles/s41467-018-063...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/SingleCell/GSE115469_GPL16791...</td>
      <td>discover-prod-datalake-v1</td>
      <td>h5ad</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626622850659</td>
      <td>Sonya A. MacParland, Ian D McGilvray</td>
      <td>2018</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```sos
sc_normal_id = list(sc_normal_query_df_filt['dataset_id'])
sc_normal_id
```




    ['GSE115469_GPL16791']



**Download single cell normal data from OmixAtlas**


```sos
url = repo_client.download_data(LiverOmixAltasID, sc_normal_id[0]).get('data')
file_name = sc_normal_id[0]+".h5ad"
os.system(f"wget -O '{file_name}' '{url}'")
```




    0



### 1.1.5 Querying Bulk RNASeq data on Liver OmixAtlas

For instance, we are required to perform RNASeq data analysis using transcriptomics data from TCGA.  
Using the curated fields we narrow down our search to find the datasets that correspond to **Liver Hepatocellular Carcinoma** in **Homo Sapiens**


```sos
liver_HCC_query = """select * FROM liveromix_atlas_files 
                    WHERE organism = 'Homo sapiens' 
                    AND disease = 'Carcinoma, Hepatocellular' """
repo_client.query_metadata(liver_HCC_query)
```

    Showing 1 - 100 of 1185 matching results





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
      <th>disease</th>
      <th>kw_drug</th>
      <th>kw_cell_line</th>
      <th>kw_cell_type</th>
      <th>publication</th>
      <th>tissue</th>
      <th>organism</th>
      <th>dataset_id</th>
      <th>kw_data_type</th>
      <th>description</th>
      <th>dataset_source</th>
      <th>curation_version</th>
      <th>platform</th>
      <th>year</th>
      <th>total_num_samples</th>
      <th>experimental_design</th>
      <th>kw_repo</th>
      <th>kw_package</th>
      <th>kw_key</th>
      <th>kw_bucket</th>
      <th>kw_filetype</th>
      <th>kw_region</th>
      <th>kw_location</th>
      <th>kw_timestamp</th>
      <th>study_id</th>
      <th>release_date</th>
      <th>kw_instrument_type</th>
      <th>publication_title</th>
      <th>kw_analysis_type</th>
      <th>is_public</th>
      <th>kw_sample_source</th>
      <th>operation</th>
      <th>data_repository</th>
      <th>kw_gene</th>
      <th>data_required</th>
      <th>file_type</th>
      <th>source_process</th>
      <th>publication_name</th>
      <th>author</th>
      <th>abstract</th>
      <th>overall_design</th>
      <th>summary</th>
      <th>type</th>
      <th>donor_age</th>
      <th>donor_sex</th>
      <th>sample_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>[Carcinoma, Hepatocellular, Hepatoblastoma]</td>
      <td>[repaglinide, vanillin, SB 216763, Nedocromil,...</td>
      <td>[Hep 3B2.1-7, Hep-G2, HLF, HuH-1, HuH-6, Huh-7...</td>
      <td>[None]</td>
      <td>https://www.biorxiv.org/content/10.1101/720243v1</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>DEPMAP_19Q4_primary-screen-mfi_liver</td>
      <td>Drug Screens</td>
      <td>Median fluorescence intensities for liver cell...</td>
      <td>DEPMAP</td>
      <td>g2</td>
      <td>Flow cytometry</td>
      <td>2019</td>
      <td>21</td>
      <td>{'categorical_variables': {'kw_curated_cell_li...</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/DrugScreens/DEPMAP_19Q4_prima...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626017983908</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>[Carcinoma, Hepatocellular, Liver Cirrhosis, N...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>https://www.ebi.ac.uk/metabolights/MTBLS105</td>
      <td>[liver, blood plasma]</td>
      <td>[Homo sapiens]</td>
      <td>MTBLS105_m_mtbls105_GC_SIM_mass_spectrometry</td>
      <td>Metabolomics</td>
      <td>This study evaluates changes in metabolite lev...</td>
      <td>Metabolights</td>
      <td>g2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>267</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Metabolomics/MTBLS105_m_mtbls...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626017990182</td>
      <td>MTBLS105</td>
      <td>2015:08:13</td>
      <td>5975C Series GC/MSD (Agilent);Agilent GC/MS, L...</td>
      <td>GC-MS based plasma metabolomics for HCC biomar...</td>
      <td>mass spectrometry</td>
      <td>true</td>
      <td>blood plasma</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>[Carcinoma, Hepatocellular, Liver Cirrhosis, F...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>https://www.ebi.ac.uk/metabolights/MTBLS17</td>
      <td>[liver, bile]</td>
      <td>[Homo sapiens]</td>
      <td>MTBLS17_m_live_mtbls17pos_metabolite profiling...</td>
      <td>Metabolomics</td>
      <td>Characterizing the metabolic changes pertainin...</td>
      <td>Metabolights</td>
      <td>g2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1050</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Metabolomics/MTBLS17_m_live_m...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626018020896</td>
      <td>MTBLS17</td>
      <td>2013:06:30</td>
      <td>UPLC Q-TOF Premier (Waters)</td>
      <td>Utilization of metabolomics to identify serum ...</td>
      <td>mass spectrometry</td>
      <td>true</td>
      <td>blood serum</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>[Hepatic Encephalopathy, Carcinoma, Hepatocell...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>https://www.ebi.ac.uk/metabolights/MTBLS582</td>
      <td>[liver]</td>
      <td>[Homo sapiens]</td>
      <td>MTBLS582_m_mtbls582_NEG_mass_spectrometry</td>
      <td>Lipidomics</td>
      <td>Obesity is tightly linked to hepatic steatosis...</td>
      <td>Metabolights</td>
      <td>g2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>12</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Metabolomics/MTBLS582_m_mtbls...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626018028414</td>
      <td>MTBLS582</td>
      <td>2018:04:26</td>
      <td>Exactive (Thermo Scientific)</td>
      <td>Global Analyses of Selective Insulin Resistanc...</td>
      <td>mass spectrometry</td>
      <td>true</td>
      <td>Hep-G2 cell</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>[Carcinoma, Hepatocellular, Liver Cirrhosis, F...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>https://www.ebi.ac.uk/metabolights/MTBLS19</td>
      <td>[liver, bile]</td>
      <td>[Homo sapiens]</td>
      <td>MTBLS19_m_neg_MTBLS19_metabolite profiling_mas...</td>
      <td>Metabolomics</td>
      <td>Although hepatocellular carcinoma (HCC) has be...</td>
      <td>Metabolights</td>
      <td>g2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>360</td>
      <td>NaN</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Metabolomics/MTBLS19_m_neg_MT...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626018026372</td>
      <td>MTBLS19</td>
      <td>2013:06:30</td>
      <td>UPLC Q-TOF Premier (Waters)</td>
      <td>LC−MS Based Serum Metabolomics for Identificat...</td>
      <td>mass spectrometry</td>
      <td>true</td>
      <td>blood serum</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>95</th>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>[BGJ-398]</td>
      <td>[Hep-G2]</td>
      <td>[None]</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>lincs_GSE70138_HEPG2_BRD-K42728290</td>
      <td>Transcriptomics</td>
      <td>Gene expression profile of Hep-G2 on treating ...</td>
      <td>LINCS</td>
      <td>g3</td>
      <td>LINCS1000</td>
      <td>2018</td>
      <td>388</td>
      <td>{'categorical_variables': {'drug': [{'name': '...</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Transcriptomics/lincs_GSE7013...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626152712356</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>true</td>
      <td>NaN</td>
      <td>{'is_normalized': 'true', 'batch_corrected_var...</td>
      <td>lincs</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>15</td>
      <td>Male</td>
      <td>Tumor</td>
    </tr>
    <tr>
      <th>96</th>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>[BKM120]</td>
      <td>[Hep-G2]</td>
      <td>[None]</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>lincs_GSE70138_HEPG2_BRD-K42191735</td>
      <td>Transcriptomics</td>
      <td>Gene expression profile of Hep-G2 on treating ...</td>
      <td>LINCS</td>
      <td>g3</td>
      <td>LINCS1000</td>
      <td>2018</td>
      <td>388</td>
      <td>{'categorical_variables': {'drug': [{'name': '...</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Transcriptomics/lincs_GSE7013...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626151670741</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>true</td>
      <td>NaN</td>
      <td>{'is_normalized': 'true', 'batch_corrected_var...</td>
      <td>lincs</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>15</td>
      <td>Male</td>
      <td>Tumor</td>
    </tr>
    <tr>
      <th>97</th>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>[Osimertinib]</td>
      <td>[Hep-G2]</td>
      <td>[None]</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>lincs_GSE70138_HEPG2_BRD-K42805893</td>
      <td>Transcriptomics</td>
      <td>Gene expression profile of Hep-G2 on treating ...</td>
      <td>LINCS</td>
      <td>g3</td>
      <td>LINCS1000</td>
      <td>2018</td>
      <td>387</td>
      <td>{'categorical_variables': {'drug': [{'name': '...</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Transcriptomics/lincs_GSE7013...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626153236520</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>true</td>
      <td>NaN</td>
      <td>{'is_normalized': 'true', 'batch_corrected_var...</td>
      <td>lincs</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>15</td>
      <td>Male</td>
      <td>Tumor</td>
    </tr>
    <tr>
      <th>98</th>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>[BX795]</td>
      <td>[Hep-G2]</td>
      <td>[None]</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>lincs_GSE70138_HEPG2_BRD-K47983010</td>
      <td>Transcriptomics</td>
      <td>Gene expression profile of Hep-G2 on treating ...</td>
      <td>LINCS</td>
      <td>g3</td>
      <td>LINCS1000</td>
      <td>2018</td>
      <td>388</td>
      <td>{'categorical_variables': {'drug': [{'name': '...</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Transcriptomics/lincs_GSE7013...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626158523464</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>true</td>
      <td>NaN</td>
      <td>{'is_normalized': 'true', 'batch_corrected_var...</td>
      <td>lincs</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>15</td>
      <td>Male</td>
      <td>Tumor</td>
    </tr>
    <tr>
      <th>99</th>
      <td>[Carcinoma, Hepatocellular]</td>
      <td>[2-((aminocarbonyl)amino)-5-(4-fluorophenyl)-3...</td>
      <td>[Hep-G2]</td>
      <td>[None]</td>
      <td>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi...</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>lincs_GSE70138_HEPG2_BRD-K51575138</td>
      <td>Transcriptomics</td>
      <td>Gene expression profile of Hep-G2 on treating ...</td>
      <td>LINCS</td>
      <td>g3</td>
      <td>LINCS1000</td>
      <td>2018</td>
      <td>388</td>
      <td>{'categorical_variables': {'drug': [{'name': '...</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/Transcriptomics/lincs_GSE7013...</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626166038957</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>true</td>
      <td>NaN</td>
      <td>{'is_normalized': 'true', 'batch_corrected_var...</td>
      <td>lincs</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>15</td>
      <td>Male</td>
      <td>Tumor</td>
    </tr>
  </tbody>
</table>
<p>100 rows × 46 columns</p>
</div>



The initial search returns Liver heaptocellular carcinoma detaasets from multiple sources and different data types.  
We are interested in liver cancer dataset from **TCGA** of type **Transcriptomics**.  

Therefore, we narrow down our search to only return datasets from source **TCGA** and type **Transcriptomics**.


```sos
liver_HCC_query = """select * FROM liveromix_atlas_files 
                    WHERE organism = 'Homo sapiens' 
                    AND disease = 'Carcinoma, Hepatocellular' 
                    AND dataset_source = 'TCGA'
                    AND kw_data_type = 'Transcriptomics'"""

liver_HCC_query_df = repo_client.query_metadata(liver_HCC_query)
liver_HCC_query_df
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
      <th>disease</th>
      <th>kw_drug</th>
      <th>kw_cell_type</th>
      <th>kw_cell_line</th>
      <th>publication</th>
      <th>tissue</th>
      <th>organism</th>
      <th>dataset_id</th>
      <th>kw_data_type</th>
      <th>description</th>
      <th>dataset_source</th>
      <th>curation_version</th>
      <th>total_num_samples</th>
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
      <td>[Carcinoma, Hepatocellular, Carcinoma, Adenoca...</td>
      <td>[Doxorubicin, sorafenib tosylate, gemcitabine,...</td>
      <td>[None]</td>
      <td>[None]</td>
      <td>https://www.cell.com/cell/fulltext/S0092-8674(...</td>
      <td>[liver]</td>
      <td>Homo sapiens</td>
      <td>LIHC_RNASeq_TCGA</td>
      <td>Transcriptomics</td>
      <td>Liver hepatocellular carcinoma RNASeq data</td>
      <td>TCGA</td>
      <td>g3</td>
      <td>424</td>
      <td>liveromix_atlas</td>
      <td>liver_atlas/data</td>
      <td>liver_atlas/data/RNAseq/LIHC_RNASeq_TCGA.gct</td>
      <td>discover-prod-datalake-v1</td>
      <td>gct</td>
      <td>us-west-2</td>
      <td>https://discover-prod-datalake-v1.s3-us-west-2...</td>
      <td>1626708777124</td>
    </tr>
  </tbody>
</table>
</div>




```sos
liver_hcc_id = list(liver_HCC_query_df['dataset_id'])
liver_hcc_id
```




    ['LIHC_RNASeq_TCGA']




```sos
url = repo_client.download_data(LiverOmixAltasID, liver_hcc_id[0]).get('data')
file_name = liver_hcc_id[0]+".gct"
os.system(f"wget -O '{file_name}' '{url}'")
```




    0



## 1.2 Outside Polly

### 1.2.1 Querying Single Cell RNASeq data on NCBI GEO

We query the datasets using the original source. The normal and tumor single cell datasets for liver patients were found on GEO.
We used the GSE Accession ID to fetch the data and metadata in the format made available by the submitters.

### Fetching liver tumor single cell data from GEO NCBI


```sos
# Fetching liver tumor data from GEO
get_data_geo('GSE125449')
```

### Fetching normal liver single cell data from GEO NCBI


```sos
# Fetching liver tumor data from GEO
get_data_geo('GSE115469')
```

### Data structure and format


```sos
cd GSE125449/ && ls -l
```

    total 45268
    -rw-r--r-- 1 polly users    26016 Jun 30 05:14 [0m[01;31mGSE125449_Set1_barcodes.tsv.gz[0m
    -rw-r--r-- 1 polly users   158399 Jun 30 05:14 [01;31mGSE125449_Set1_genes.tsv.gz[0m
    -rw-r--r-- 1 polly users 26707779 Jun 30 05:14 [01;31mGSE125449_Set1_matrix.mtx.gz[0m
    -rw-r--r-- 1 polly users    33526 Jun 30 05:14 [01;31mGSE125449_Set1_samples.txt.gz[0m
    -rw-r--r-- 1 polly users    23861 Jun 30 05:14 [01;31mGSE125449_Set2_barcodes.tsv.gz[0m
    -rw-r--r-- 1 polly users   153470 Jun 30 05:14 [01;31mGSE125449_Set2_genes.tsv.gz[0m
    -rw-r--r-- 1 polly users 19204961 Jun 30 05:14 [01;31mGSE125449_Set2_matrix.mtx.gz[0m
    -rw-r--r-- 1 polly users    31132 Jun 30 05:14 [01;31mGSE125449_Set2_samples.txt.gz[0m



```sos
cd ../GSE115469/ && ls -l 
```

    total 105536
    -rw-r--r-- 1 polly users 108065293 Jun 30 05:15 [0m[01;31mGSE115469_Data.csv.gz[0m


### 1.2.2 Querying Bulk RNASeq data on TCGA

### Fetching mRNA data from TCGA

```R
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM",
                  legacy = FALSE)
GDCdownload(query)
data <- GDCprepare(query)
```

**Time taken to fetch data: ~ 10 mins**

```
Downloading data for project TCGA-LIHC
GDCdownload will download 424 files. A total of 207.74981 MB
Downloading as: Tue_Jun_29_11_38_54_2021.tar.gz
Downloading: 210 MB
|=====================================================================|100%                      
Completed after 3 m 
Starting to add information to samples
 => Add clinical information to samples
 => Adding TCGA molecular information from marker papers
 => Information will have prefix 'paper_' 
From the 60483 genes we couldn't map 3881
user  system elapsed 
265.095  61.711 519.107
```

### Data Structure and format

```
class: RangedSummarizedExperiment 
dim: 56602 424 
metadata(1): data_release
assays(1): HTSeq - FPKM
rownames(56602): ENSG00000000003 ENSG00000000005 ... ENSG00000281912
  ENSG00000281920
rowData names(3): ensembl_gene_id external_gene_name original_ensembl_gene_id
colnames(424): TCGA-DD-AAE4-01A-11R-A41C-07 TCGA-CC-A3M9-01A-11R-A213-07 ...
  TCGA-DD-AAE0-01A-11R-A41C-07 TCGA-DD-A3A8-11A-11R-A22L-07
colData names(64): barcode patient ... released sample.aux
```

# 2. Single Cell data integration

### 2.1 Reading single cell datasets

Reading h5ad as an anndata object and converting it to a Seurat object contained assay, gene and sample metadata


```sos
suppressMessages(library("Seurat"))
suppressMessages(library("anndata"))
suppressMessages(library("dplyr"))
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("clustermole"))
options(warn=-1)
```


```sos
#Create a Seurat object 
get_seurat_obj <- function(adata) {
  
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
  
  return(seurat)
}
```


```sos
adata_normal <- read_h5ad('GSE115469_GPL16791.h5ad')
adata_tumor_1 <- read_h5ad('GSE125449_GPL18573.h5ad')
adata_tumor_2 <- read_h5ad('GSE125449_GPL20301.h5ad')
```


```sos
normal_seurat <- get_seurat_obj(adata_normal)
tumor_seurat_1 <- get_seurat_obj(adata_tumor_1)
tumor_seurat_2 <- get_seurat_obj(adata_tumor_2)
```


```sos
#Add condition column to metadata
normal_seurat@meta.data$condition <- c('Normal')
tumor_seurat_1@meta.data$condition <- c('Tumor')
tumor_seurat_2@meta.data$condition <- c('Tumor')

normal_seurat <- PercentageFeatureSet(normal_seurat, pattern = "^MT-", col.name = "percent.mt")
tumor_seurat_1 <- PercentageFeatureSet(tumor_seurat_1, pattern = "^MT-", col.name = "percent.mt")
tumor_seurat_2 <- PercentageFeatureSet(tumor_seurat_2, pattern = "^MT-", col.name = "percent.mt")
```

### 2.2 Identification of cell types for normal single cell dataset


```sos
# find markers for every cluster compared to all remaining cells, report only the positive ones
Idents(object = normal_seurat) <- normal_seurat@meta.data$clusters
liver.markers <- FindAllMarkers(normal_seurat,min.pct = 0.25, logfc.threshold = 0.25)
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
```


```sos
cell_type_clust <- aggregate(clusters ~ Type, clust_map, paste, collapse=',')
cell_type_clust
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
#Add cell types to seurat metadata
row.meta <- rownames(normal_seurat@meta.data)
normal_seurat@meta.data <- merge(normal_seurat@meta.data, clust_map, by = 'clusters')
rownames(normal_seurat@meta.data) <- row.meta
```

### 2.3 Single cell data integration 


```sos
seurat <- merge(normal_seurat, y = c(tumor_seurat_1, tumor_seurat_2))

seurat.list <- SplitObject(seurat, split.by = "condition")
```


```sos
#Clear memory 
rm(normal_seurat, tumor_seurat_1, tumor_seurat_2)
```

Performing normalisation and scaling using `SCTransform`


```sos
seurat.list[["Normal"]] <- SCTransform(seurat.list[["Normal"]], vars.to.regress = 'percent.mt')
seurat.list[["Tumor"]] <- SCTransform(seurat.list[["Tumor"]], vars.to.regress = 'percent.mt')
```


```sos
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
liver.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = features)
```


```sos
liver.combined.sct <- IntegrateData(anchorset = liver.anchors, normalization.method = "SCT")
```

    Merging dataset 1 into 2
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    


Compute PCA and UMAP using first 40 principal components


```sos
liver.combined.sct <- RunPCA(liver.combined.sct, verbose = FALSE)
liver.combined.sct <- RunUMAP(liver.combined.sct, reduction = "pca", dims = 1:40)
```

Perform SNN clustering of cells using first 40 principal components


```sos
liver.combined.sct <- FindNeighbors(liver.combined.sct, dims = 1:40)
liver.combined.sct <- FindClusters(liver.combined.sct, resolution = 0.8)
```

Plot the UMAP clustering of cell type and cluster populations


```sos
DimPlot(liver.combined.sct, reduction = "umap", group.by = 'seurat_clusters', cols = DiscretePalette(30), label = TRUE)
DimPlot(liver.combined.sct, reduction = "umap", group.by = 'Type', cols = DiscretePalette(19))
```


Marker analysis of fibroblast population


```sos
# find markers for fibroblast clusters compared to all remaining cells, report only the positive ones
DefaultAssay(liver.combined.sct) <- "RNA"
Idents(object = liver.combined.sct) <- liver.combined.sct@meta.data$seurat_clusters
clus7.markers <- FindMarkers(liver.combined.sct, ident.1 = '7', min.pct = 0.25)
clus15markers <- FindMarkers(liver.combined.sct, ident.1 = '15', min.pct = 0.25)
```

Selecting markers from Fibroblast sub populations


```sos
VlnPlot(liver.combined.sct, features = "COL1A1", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "COL1A2", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "COL3A1", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "BGN", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "LUM", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "SPARCL1", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "GJA4", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "OAZ2", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "ADAMTS4", idents = c("7","15"))
VlnPlot(liver.combined.sct, features = "GPR4", idents = c("7","15"))
```

We observed that fibrogenic genes (e.g., COL1A1, COL1A2, COL3A1 ) and proteoglycan encoding genes (e.g., LUM, BGN ) were overexpressed in cluster 15. Cluster 7, however showed high overall expression of SPARCL1, GJA4. 

In order to assess the clinical significance bulk RNASeq data was fetched from Liver OmixAtlas for Heaptocellular Carcinoma.

# 3. Survival analysis using bulk RNASeq data

### Reading the mRNA data 


```sos
# Parsing mRNA expression data 
gct = parse('./LIHC_RNASeq_TCGA.gct')

# column metadata
clinical_data = gct.col_metadata_df

# expression data
expression_data = gct.data_df
```

### Getting paired normal and tumor samples


```sos
clinical_data_paired = clinical_data[clinical_data.duplicated(subset = 'patient', keep = False)]
clinical_data_paired.head()
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
      <th>chd</th>
      <th>patient</th>
      <th>barcode</th>
      <th>sample</th>
      <th>shortLetterCode</th>
      <th>definition</th>
      <th>sample_submitter_id</th>
      <th>sample_type_id</th>
      <th>sample_id</th>
      <th>sample_type</th>
      <th>state</th>
      <th>initial_weight</th>
      <th>pathology_report_uuid</th>
      <th>submitter_id</th>
      <th>oct_embedded</th>
      <th>is_ffpe</th>
      <th>tissue_type</th>
      <th>synchronous_malignancy</th>
      <th>ajcc_pathologic_stage</th>
      <th>tumor_stage</th>
      <th>last_known_disease_status</th>
      <th>tissue_or_organ_of_origin</th>
      <th>primary_diagnosis</th>
      <th>prior_malignancy</th>
      <th>prior_treatment</th>
      <th>ajcc_staging_system_edition</th>
      <th>ajcc_pathologic_t</th>
      <th>morphology</th>
      <th>ajcc_pathologic_n</th>
      <th>ajcc_pathologic_m</th>
      <th>classification_of_tumor</th>
      <th>diagnosis_id</th>
      <th>icd_10_code</th>
      <th>site_of_resection_or_biopsy</th>
      <th>tumor_grade</th>
      <th>progression_or_recurrence</th>
      <th>alcohol_history</th>
      <th>exposure_id</th>
      <th>weight</th>
      <th>height</th>
      <th>bmi</th>
      <th>...</th>
      <th>therapy_regimen.drug_4</th>
      <th>therapy_regimen_other.drug_4</th>
      <th>total_dose.drug_4</th>
      <th>tx_on_clinical_trial.drug_4</th>
      <th>AFP&gt;300</th>
      <th>Baylor GCC whole genome seq</th>
      <th>bcr_drug_barcode.drug_5</th>
      <th>bcr_drug_uuid.drug_5</th>
      <th>form_completion_date.drug_5</th>
      <th>pharmaceutical_therapy_drug_name.drug_5</th>
      <th>clinical_trial_drug_classification.drug_5</th>
      <th>pharmaceutical_therapy_type.drug_5</th>
      <th>pharmaceutical_tx_started_days_to.drug_5</th>
      <th>pharmaceutical_tx_ongoing_indicator.drug_5</th>
      <th>pharmaceutical_tx_ended_days_to.drug_5</th>
      <th>treatment_best_response.drug_5</th>
      <th>days_to_stem_cell_transplantation.drug_5</th>
      <th>pharm_regimen.drug_5</th>
      <th>pharm_regimen_other.drug_5</th>
      <th>pharma_adjuvant_cycles_count.drug_5</th>
      <th>pharma_type_other.drug_5</th>
      <th>pharmaceutical_tx_dose_units.drug_5</th>
      <th>pharmaceutical_tx_total_dose_units.drug_5</th>
      <th>prescribed_dose.drug_5</th>
      <th>regimen_number.drug_5</th>
      <th>route_of_administration.drug_5</th>
      <th>stem_cell_transplantation.drug_5</th>
      <th>stem_cell_transplantation_type.drug_5</th>
      <th>therapy_regimen.drug_5</th>
      <th>therapy_regimen_other.drug_5</th>
      <th>total_dose.drug_5</th>
      <th>tx_on_clinical_trial.drug_5</th>
      <th>kw_curated_disease</th>
      <th>kw_curated_drug</th>
      <th>kw_curated_tissue</th>
      <th>kw_curated_cell_type</th>
      <th>kw_curated_cell_line</th>
      <th>kw_curated_genetic_mod_type</th>
      <th>kw_curated_modified_gene</th>
      <th>kw_curated_gene</th>
    </tr>
    <tr>
      <th>cid</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-BC-A10Q-01A-11R-A131-07</th>
      <td>TCGA-BC-A10Q</td>
      <td>TCGA-BC-A10Q-01A-11R-A131-07</td>
      <td>TCGA-BC-A10Q-01A</td>
      <td>TP</td>
      <td>Primary solid Tumor</td>
      <td>TCGA-BC-A10Q-01A</td>
      <td>1</td>
      <td>bf225dc1-4309-4304-849e-cbcc23c8442c</td>
      <td>Primary Tumor</td>
      <td>released</td>
      <td>340</td>
      <td>10117C34-85C2-49DE-B2AD-F8D0D4570DA9</td>
      <td>TCGA-BC-A10Q</td>
      <td>false</td>
      <td>FALSE</td>
      <td>Not Reported</td>
      <td>No</td>
      <td>NaN</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Liver</td>
      <td>Hepatocellular carcinoma, NOS</td>
      <td>no</td>
      <td>No</td>
      <td>5th</td>
      <td>T2</td>
      <td>8170/3</td>
      <td>NX</td>
      <td>MX</td>
      <td>not reported</td>
      <td>c896279a-640a-57f8-91c5-4d1b674cdcca</td>
      <td>C22.0</td>
      <td>Liver</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Not Reported</td>
      <td>41801e4c-8d70-5266-82ae-382fc1dd4a20</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>PROGRESSION</td>
      <td>[Not Applicable]</td>
      <td>[Not Available]</td>
      <td>[Not Available]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Carcinoma, Hepatocellular</td>
      <td>Doxorubicin</td>
      <td>liver</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
    </tr>
    <tr>
      <th>TCGA-BC-A10Q-11A-11R-A131-07</th>
      <td>TCGA-BC-A10Q</td>
      <td>TCGA-BC-A10Q-11A-11R-A131-07</td>
      <td>TCGA-BC-A10Q-11A</td>
      <td>NT</td>
      <td>Solid Tissue Normal</td>
      <td>TCGA-BC-A10Q-11A</td>
      <td>11</td>
      <td>5d3bf988-3331-4412-bc17-b8f4650a8623</td>
      <td>Solid Tissue Normal</td>
      <td>released</td>
      <td>340</td>
      <td>NaN</td>
      <td>TCGA-BC-A10Q</td>
      <td>false</td>
      <td>FALSE</td>
      <td>Not Reported</td>
      <td>No</td>
      <td>NaN</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Liver</td>
      <td>Hepatocellular carcinoma, NOS</td>
      <td>no</td>
      <td>No</td>
      <td>5th</td>
      <td>T2</td>
      <td>8170/3</td>
      <td>NX</td>
      <td>MX</td>
      <td>not reported</td>
      <td>c896279a-640a-57f8-91c5-4d1b674cdcca</td>
      <td>C22.0</td>
      <td>Liver</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Not Reported</td>
      <td>41801e4c-8d70-5266-82ae-382fc1dd4a20</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>PROGRESSION</td>
      <td>[Not Applicable]</td>
      <td>[Not Available]</td>
      <td>[Not Available]</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Normal</td>
      <td>Doxorubicin</td>
      <td>liver</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
    </tr>
    <tr>
      <th>TCGA-BC-A10R-01A-11R-A131-07</th>
      <td>TCGA-BC-A10R</td>
      <td>TCGA-BC-A10R-01A-11R-A131-07</td>
      <td>TCGA-BC-A10R-01A</td>
      <td>TP</td>
      <td>Primary solid Tumor</td>
      <td>TCGA-BC-A10R-01A</td>
      <td>1</td>
      <td>995180f2-5aa7-47b9-bb3f-585f94b2457a</td>
      <td>Primary Tumor</td>
      <td>released</td>
      <td>180</td>
      <td>65BE2C95-CADA-44FC-9313-14529AFF773E</td>
      <td>TCGA-BC-A10R</td>
      <td>false</td>
      <td>FALSE</td>
      <td>Not Reported</td>
      <td>No</td>
      <td>NaN</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Liver</td>
      <td>Hepatocellular carcinoma, NOS</td>
      <td>no</td>
      <td>No</td>
      <td>5th</td>
      <td>T3</td>
      <td>8170/3</td>
      <td>NX</td>
      <td>MX</td>
      <td>not reported</td>
      <td>ef53ae43-7ab9-546f-b40c-71018f4bf37a</td>
      <td>C22.0</td>
      <td>Liver</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Not Reported</td>
      <td>d17d5f9e-bcd9-5955-baaa-f03189d6006a</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>No</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Carcinoma, Hepatocellular</td>
      <td>none</td>
      <td>liver</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
    </tr>
    <tr>
      <th>TCGA-BC-A10R-11A-11R-A131-07</th>
      <td>TCGA-BC-A10R</td>
      <td>TCGA-BC-A10R-11A-11R-A131-07</td>
      <td>TCGA-BC-A10R-11A</td>
      <td>NT</td>
      <td>Solid Tissue Normal</td>
      <td>TCGA-BC-A10R-11A</td>
      <td>11</td>
      <td>1ea00561-3d1d-49d0-a3f7-d34eb8fec235</td>
      <td>Solid Tissue Normal</td>
      <td>released</td>
      <td>1620</td>
      <td>NaN</td>
      <td>TCGA-BC-A10R</td>
      <td>false</td>
      <td>FALSE</td>
      <td>Not Reported</td>
      <td>No</td>
      <td>NaN</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Liver</td>
      <td>Hepatocellular carcinoma, NOS</td>
      <td>no</td>
      <td>No</td>
      <td>5th</td>
      <td>T3</td>
      <td>8170/3</td>
      <td>NX</td>
      <td>MX</td>
      <td>not reported</td>
      <td>ef53ae43-7ab9-546f-b40c-71018f4bf37a</td>
      <td>C22.0</td>
      <td>Liver</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Not Reported</td>
      <td>d17d5f9e-bcd9-5955-baaa-f03189d6006a</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>No</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Normal</td>
      <td>none</td>
      <td>liver</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
    </tr>
    <tr>
      <th>TCGA-BC-A10T-01A-11R-A131-07</th>
      <td>TCGA-BC-A10T</td>
      <td>TCGA-BC-A10T-01A-11R-A131-07</td>
      <td>TCGA-BC-A10T-01A</td>
      <td>TP</td>
      <td>Primary solid Tumor</td>
      <td>TCGA-BC-A10T-01A</td>
      <td>1</td>
      <td>d2f6cea0-ecec-49f1-a851-e332e79ce098</td>
      <td>Primary Tumor</td>
      <td>released</td>
      <td>1000</td>
      <td>D8103C33-F303-4D97-81AF-CC3D06E79285</td>
      <td>TCGA-BC-A10T</td>
      <td>false</td>
      <td>FALSE</td>
      <td>Not Reported</td>
      <td>No</td>
      <td>NaN</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Liver</td>
      <td>Hepatocellular carcinoma, NOS</td>
      <td>no</td>
      <td>No</td>
      <td>5th</td>
      <td>T4</td>
      <td>8170/3</td>
      <td>NX</td>
      <td>MX</td>
      <td>not reported</td>
      <td>b5925a40-cb9f-5954-9cfb-c7ab86404e60</td>
      <td>C22.0</td>
      <td>Liver</td>
      <td>not reported</td>
      <td>not reported</td>
      <td>Not Reported</td>
      <td>a1d6658a-accb-5864-a6b2-31e176486998</td>
      <td>86.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Carcinoma, Hepatocellular</td>
      <td>none</td>
      <td>liver</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
      <td>none</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 301 columns</p>
</div>




```sos
#Subset the expression data using the paired normal and tumor samples
expression_data_paired = expression_data[expression_data.columns.intersection(clinical_data_paired.index)]
```


```sos
expression_data_paired.to_csv('TCGA_LIHC_paired.csv', sep = '\t')
clinical_data_paired.to_csv('TCGA_LIHC_clinical.csv', sep = '\t')
```

### Preparing TCGA clinical data for survival analysis


```sos
#Load R packages
suppressMessages(library("readr"))
suppressMessages(library("gplots"))
suppressMessages(library("survival"))
suppressMessages(library("survminer"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gProfileR"))
suppressMessages(library("genefilter"))
```


```sos
expression_data <- as.data.frame(read_delim("TCGA_LIHC_paired.csv", "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols()))
rownames(expression_data) <- expression_data[,1]
expression_data <- expression_data[,-1]

clinical_data <- as.data.frame(read_delim("TCGA_LIHC_clinical.csv", "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols()))
rownames(clinical_data) <- clinical_data[,1]
clinical_data <- clinical_data[,-1]
```


```sos
# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clinical_data$deceased = clinical_data$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_data$overall_survival = ifelse(clinical_data$deceased,
                                   clinical_data$days_to_death,
                                   clinical_data$days_to_last_follow_up)

# show first 10 samples
clinical_data_survival = clinical_data[,c('vital_status','deceased','days_to_death','days_to_last_follow_up','overall_survival')]
head(clinical_data_survival)
```


<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>vital_status</th><th scope=col>deceased</th><th scope=col>days_to_death</th><th scope=col>days_to_last_follow_up</th><th scope=col>overall_survival</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TCGA-BC-A10Q-01A-11R-A131-07</th><td>Dead</td><td>TRUE</td><td>1135</td><td>NA</td><td>1135</td></tr>
	<tr><th scope=row>TCGA-BC-A10Q-11A-11R-A131-07</th><td>Dead</td><td>TRUE</td><td>1135</td><td>NA</td><td>1135</td></tr>
	<tr><th scope=row>TCGA-BC-A10R-01A-11R-A131-07</th><td>Dead</td><td>TRUE</td><td> 308</td><td>NA</td><td> 308</td></tr>
	<tr><th scope=row>TCGA-BC-A10R-11A-11R-A131-07</th><td>Dead</td><td>TRUE</td><td> 308</td><td>NA</td><td> 308</td></tr>
	<tr><th scope=row>TCGA-BC-A10T-01A-11R-A131-07</th><td>Dead</td><td>TRUE</td><td> 837</td><td>NA</td><td> 837</td></tr>
	<tr><th scope=row>TCGA-BC-A10T-11A-11R-A131-07</th><td>Dead</td><td>TRUE</td><td> 837</td><td>NA</td><td> 837</td></tr>
</tbody>
</table>



### Examining the association of cluster specific gene sets with clinical data

1. Specific Cluster 24 gene set G24 : ( COL1A, CFH, LUM and PDGFRA ). 
2. Specific Cluster 8 gene set G3 : ( SPARCL1, GJA4, ADAMTS4 and GPR4 ). 


```sos
survival_gene_exp = function(gene_name) {
    
    expression_mat = t(expression_data)
    # get the expression values for the selected gene
    clinical_data$gene_name = expression_mat[rownames(clinical_data), gene_name]

    # find the median value of the gene and print it
    median_value = median(clinical_data$gene_name)
    
    # divide patients in two groups, up and down regulated.
    # if the patient expression is greater or equal to them median we put it
    # among the "up-regulated", otherwise among the "down-regulated"
    clinical_data$gene = ifelse(clinical_data$gene_name >= median_value, "UP", "DOWN")

    # we can fit a survival model, like we did in the previous section
    fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clinical_data)
    
    return(list(fit, clinical_data))
    
}
```

### Survival Analysis of Cluster 8 genes

We observe that overexpression of SPARCL1 gene showed a favourable prognosis in liver cancer patients. 


```sos
# and finally, we produce a Kaplan-Meier plot
surv = survival_gene_exp("SPARCL1")
ggsurvplot(surv[[1]], data=surv[[2]], pval=T, risk.table=T, title=paste("SPARCL1"))

surv = survival_gene_exp("GJA4")
ggsurvplot(surv[[1]], data=surv[[2]], pval=T, risk.table=T, title=paste("GJA4"))

surv = survival_gene_exp("ADAMTS4")
ggsurvplot(surv[[1]], data=surv[[2]], pval=T, risk.table=T, title=paste("ADAMTS4"))

surv = survival_gene_exp("GPR4")
ggsurvplot(surv[[1]], data=surv[[2]], pval=T, risk.table=T, title=paste("GPR4"))
```

## Appendix

### Package installations

*NOTE : For R packages, the cells must be run in the following order as shown below*


```sos
!pip3 install cmapPy --user
!pip3 install GEOparse --user
!pip3 install sh --user
!pip3 install wget --user
```


```sos
!sudo R -e 'BiocManager::install(c("genefilter", "SingleCellExperiment","GSEABase", "GSVA", "singscore"))'
```


```sos
!sudo R -e 'install.packages(c("anndata","clustermole","gridExtra"), repos = "https://cloud.r-project.org/")'
```


```sos
!sudo R -e 'install.packages(c("survival", "survminer", "gProfileR", "RColorBrewer", "gplots"), repos = "https://cloud.r-project.org/")'
```
