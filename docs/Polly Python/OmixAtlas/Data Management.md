### save_to_workspace()
Save datasets from an OmixAtlas to a workspace

Example to save the dataset_id 'GSE101127_GPL1355' from repo_id 1615965444377 to a workspace_id 8025 in a folder named 'data'
```
omixatlas.save_to_workspace('1615965444377', 'GSE101127_GPL1355', 8025, 'data')
```

### add_datasets()

Data-admin with appropriate repository level access can ingest data to an OmixAtlas using the following function:- 
```
add_datasets(repo_id (int/str), source_folder_path (dict), destination_folder (str) (optional), priority (str) (optiona))
```

input: 

`repo_id`: This is the repository ID to which ingestion should be done

`source_folder_path`: This is the dictionary with keys data and metadata. The corresponding value pairs should be the folder containing the file (gct, h5ad, vcf, mmcif etc) for data and folder containing json of dataset level metadata for metadata.

`destination_folder` (optional): This is the folder within S3 when data gets pushed

`priority` (optional): This is the priority flag as per ingestion is being done. Default is 'medium'

output: 

Status of file upload for each dataset in a dataframe

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()

repo_id = "1657110718820"
source_folder_path_data = "/import/data_final"
source_folder_path_metadata = "/import/metadata_final"
destination_folder = "220707-1426"
priority = "high"
source_folder_path = {"data":source_folder_path_data, "metadata":source_folder_path_metadata}
omixatlas.add_datasets(repo_id, source_folder_path, destination_folder, priority)
```


One of the ingestion case studies where the user wants to ingest data available in an existing OmixAtlas to a new OmixAtlas is described in this [notebook](https://github.com/ElucidataInc/polly-python/blob/main/Ingest/Data_Ingestion_CaseID_1.ipynb)

### delete_datasets()
Data-admin with appropriate repository level access can delete data from an OmixAtlas using the following function:- 

```
delete_datasets(repo_id: int/str, dataset_ids: list<string>):
```

input:

`repo_id`: (int/str) This is the repository ID from which dataset should be deleted

`dataset_ids`: (list<string>) dataset_ids that users want to delete

output:

Status of file deletion for each dataset in a dataframe

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()

repo_id = "1657110718820"
dataset_ids = ["GSE12345_GPL123", "GSE56789_GPL456"]

omixatlas.delete_datasets(repo_id, dataset_ids)
```

### dataset_metadata_template()
In order to ingest the dataset level metadata appropriately in the OmixAtlas, the user needs to ensure the json files contains the keys as per the dataset level schema of the OmixAtlas. A template of the keys can be fetched using the following function:-

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()

dataset_metadata_template(repo_id: int/str)
```
input:

repo_id: This is the repository ID for which user wants to generate template of dataset level metadata.

output: 

A dictionary with all the field names which should be present in the dataset level metadata

### update_datasets()
```update_datasets(repo_id (int/str), source_folder_path (dict), destination_folder (str) (optional), priority (str) (optiona))```

input: 

```repo_id```: This is the repository ID to which ingestion should be done

```source_folder_path```: This is the dictionary with keys ```data``` and ```metadata```. The corresponding value pairs should be the folder containing the file (gct, h5ad, vcf, mmcif etc) for data and folder containing json of dataset level metadata for metadata. The user can provide either data or metadata or both depending on which files the user wants to update. 

```destination_folder``` (optional): This is the folder within S3 when data gets pushed

```priority``` (optional): This is the priority flag as per ingestion is being done. Default is 'medium'

output: 

Status of file update for each dataset in a dataframe

```
#incase the user wants to update both metadata and data
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()

repo_id = "1657110718820"
source_folder_path_data = "/import/data_final"
source_folder_path_metadata = "/import/metadata_final"
destination_folder = "220707-1426"
priority = "high"
source_folder_path = {"data":source_folder_path_data, "metadata":source_folder_path_metadata}
omixatlas.update_datasets(repo_id, source_folder_path, destination_folder, priority)

#incase the user wants to update only the metadata
source_folder_path = {"metadata":source_folder_path_metadata}
omixatlas.update_datasets(repo_id, source_folder_path, destination_folder, priority)

#incase the user wants to update only data 
source_folder_path = {"data":source_folder_path_data}
omixatlas.update_datasets(repo_id, source_folder_path, destination_folder, priority)
```

Note: 
1. metadata/data files will not be updated if they have never been uploaded onto the omixatlas before. 
