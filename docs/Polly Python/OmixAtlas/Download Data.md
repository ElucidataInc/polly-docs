::: polly.Download
    options:
      show_source: false

## Examples

### download_data()

#### Downloading .gct and opening it in a data frame

```py
dataset_id = "GSE100003_GPL15207" #dataset which user wants to download.
repo_key = 9 OR "geo" #repo_id OR repo_name from which dataset should be downloaded from.
file_name = f"{dataset_id}.gct"
data = client.download_data(repo_key, dataset_id)
url = data.get('data').get('attributes').get('download_url')
status = os.system(f"wget -O '{file_name}' '{url}'")
if status == 0:
    print("Downloaded data successfully")
else:
    raise Exception("Download not successful")
```

```py
# In order to parse the .gct data, a python package called cmapPy can be used in the following manner.
import pandas as pd
import cmapPy
from cmapPy.pandasGEXpress.parse_gct import parse

gct_obj = parse(file_name) # Parse the file to create a gct object
df_real = gct_obj.data_df # Extract the dataframe from the gct object
col_metadata = gct_obj.col_metadata_df # Extract the column metadata from the gct object
row_metadata = gct_obj.row_metadata_df # Extract the row metadata from the gct object
```

#### Downloading .h5ad file and opening it in a data frame

```py
dataset_id = "GSE121001_GPL19057" #dataset which user wants to download.
repo_key = 17 OR "sc_data_lake" #repo_id OR repo_name from which dataset should be downloaded from.
file_name = f"{dataset_id}.h5ad"
data = client.download_data(repo_key, dataset_id)
url = data.get('data').get('attributes').get('download_url')
status = os.system(f"wget -O '{file_name}' '{url}'")
if status == 0:
    print("Downloaded data successfully")
else:
    raise Exception("Download not successful")
```

```py
# In order to parse the .h5ad data, a python package called scanpy can be used in the following manner.
import pandas as pd
import scanpy
data = sc.read_h5ad(file_name)
obs = data.obs.head()
var = data.var.head()
```

#### Downloading vcf files

```py
dataset_id = "gnomad_v2.1.1_genome_TP53" #dataset which user wants to download.
repo_key = 1628836648493 OR "gnomad" #repo_id OR repo_name from which dataset should be downloaded from.
file_name = f"{dataset_id}.vcf"
data = client.download_data(repo_key, dataset_id)
url = data.get('data').get('attributes').get('download_url')
status = os.system(f"wget -O '{file_name}' '{url}'")
if status == 0:
    print("Downloaded data successfully")
else:
    raise Exception("Download not successful")</code></pre>
```

The downloaded vcf file can be further analysed using the docker environment containing Hail package on Polly.

### download_metadata