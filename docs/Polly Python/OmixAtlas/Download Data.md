### download_data()
To download any dataset, the following function can be used to get the signed URL of the dataset.

<pre><code>omixatlas.download_data("repo_key", "[dataset_id]")</code></pre>

`repo_key` (str) repo_id OR repo_name from where the data needs to be downloaded.

`dataset_id` (str) dataset_id which the user wants to download.

The <code>[repo_name OR repo_id]</code> of an OmixAtlas can be identified by calling the <code>get_all_omixatlas()</code> function. The <code>[dataset_id]</code> can be obtained by querying the metadata at the dataset level using <code>query_metadata("[query written in SQL]")</code>.

The output of this function is a *signed URL*. The data can be downloaded by clicking on this URL.

> **_NOTE:_** This signed URL expires after 60 minutes from when it is generated.

<br>The output data is in gct/h5ad/vcf/mmcif format. This data can be parsed into a data frame for better accessibility using the following code:

#### Downloading .gct and opening it in a data frame
<pre><code>dataset_id = "GSE100003_GPL15207" #dataset which user wants to download.
repo_key = 9 OR "geo" #repo_id OR repo_name from which dataset should be downloaded from.
file_name = f"{dataset_id}.gct"
data = client.download_data(repo_key, dataset_id)
url = data.get('data').get('attributes').get('download_url')
status = os.system(f"wget -O '{file_name}' '{url}'")
if status == 0:
    print("Downloaded data successfully")
else:
    raise Exception("Download not successful")
</code></pre>

In order to parse the .gct data, a python package called cmapPy can be used in the following manner.

<pre><code>import pandas as pd
import cmapPy
from cmapPy.pandasGEXpress.parse_gct import parse

gct_obj = parse(file_name) # Parse the file to create a gct object
df_real = gct_obj.data_df # Extract the dataframe from the gct object
col_metadata = gct_obj.col_metadata_df # Extract the column metadata from the gct object
row_metadata = gct_obj.row_metadata_df # Extract the row metadata from the gct object
</code></pre>

#### Downloading .h5ad file and opening it in a data frame
<pre><code>dataset_id = "GSE121001_GPL19057" #dataset which user wants to download.
repo_key = 17 OR "sc_data_lake" #repo_id OR repo_name from which dataset should be downloaded from.
file_name = f"{dataset_id}.h5ad"
data = client.download_data(repo_key, dataset_id)
url = data.get('data').get('attributes').get('download_url')
status = os.system(f"wget -O '{file_name}' '{url}'")
if status == 0:
    print("Downloaded data successfully")
else:
    raise Exception("Download not successful")
</code></pre>

In order to parse the .h5ad data, a python package called scanpy can be used in the following manner.

<pre><code>import pandas as pd
import scanpy
data = sc.read_h5ad(file_name)
obs = data.obs.head()
var = data.var.head()
</code></pre>

In order to get started with analysis of single cell data on Polly, users can refer to this [notebook](https://github.com/ElucidataInc/polly-python/blob/main/consumption_starter_notebooks/SingleCell-polly-python.ipynb) hosted on our github.

#### Downloading vcf files
<pre><code>dataset_id = "gnomad_v2.1.1_genome_TP53" #dataset which user wants to download.
repo_key = 1628836648493 OR "gnomad" #repo_id OR repo_name from which dataset should be downloaded from.
file_name = f"{dataset_id}.vcf"
data = client.download_data(repo_key, dataset_id)
url = data.get('data').get('attributes').get('download_url')
status = os.system(f"wget -O '{file_name}' '{url}'")
if status == 0:
    print("Downloaded data successfully")
else:
    raise Exception("Download not successful")</code></pre>

The downloaded vcf file can be further analysed using the docker environment containing Hail package on Polly.

### download_metadata()

This function is used to download the dataset level metadata into a json file. The keys of the json file is kept as original_name in the schema.
```
1 from polly.omixatlas import OmixAtlas
2 omixatlas = OmixAtlas (AUTH_TOKEN)
3 omixatlas.download_metadata(repo_key, dataset_id, file_path)
``` 
Argument description:-

repo_key (str): repo_name/repo_id of the repository where data exists.
dataset_id (str): dataset_id of the dataset for which the metadata to be downloaded
file_path (str): the system destination path where the dataset level metadata should be downloaded.
