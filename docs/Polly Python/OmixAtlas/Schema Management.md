To enable users to interact with the schema of a particular OmixAtlas, functions for visualizing, updating and inserting schema is released. Updating and inseting schema is allowed for users who have data-admin credentials only.

### get_schema()
Use `get_schema(repo_key, schema_level: list (optional), source: str (optional), data_type: str (optional), return_type: (str) (optional)` to extract the schema of an OmixAtlas.

`repo_key` repo_id OR repo_name. This is a mandatory field. 

`schema_level` (optional) Table names for the OmixAtlas as per response of `query_metadata` function for the following query: `SHOW TABLES IN <repo_name>` For example:

`datasets`, `samples` and `features` for gct files in geo OmixAtlas (`geo`)

`samples_singlecell` and `features_singlecell` for h5ad files in Single Cell OmixAtlas (`sc_data_lake`)

`variant_data` for vcf files in gnomad OmixAtlas (`gnomad`)

This is an optional parameter, and by default schema for all the tables will be fetched.

`source` (optional) if source specific schema is ingested in the OA, then this field can be used to extract that. 

`data_type` (optional) if data type specific schema is ingested in the OA, then this field can be used to extract that 

`return_type` (optional) takes two inputs dataframe and dict with the default value to be dataframe

Example to fetch dataset, sample and feature level schema for all datatypes from all sources in GEO Omixatlas in dictionary format

```
schema = omixatlas.get_schema("geo", ['datasets', 'samples', 'features'], "all", "all", "dict")

to fetch the dictionary with entire payload of dataset level metadata,

`schema.datasets`

to fetch the dictionary with entire payload of sample level metadata,

`schema.samples`

to fetch the dictionary with entire payload of feature level metadata,

`schema.features`
```

Similarly, get_schema will give dataframe output for the following:-

```
schema = omixatlas.get_schema("geo", ['datasets', 'samples', 'features'], "all", "all", "dataframe")

to fetch the dataframe with summary of dataset level metadata,

`schema.datasets`

to fetch the dataframe with summary of sample level metadata,

`schema.samples`

to fetch the dataframe with entire payload of feature level metadata,

`schema.features`
```

### update_schema()
Use `update_schema(repo_key, payload)` to update the existing schema of an OmixAtlas.

```
omixatlas.update_schema(repo_key, payload)
```

`repo_key`: (str) repo_id OR repo_name. This is a mandatory field.

`payload`: (dict) The payload is a JSON file which should be as per the structure defined for schema. Only data-admin will have the authentication to update the schema.

`payload` can be loaded from the JSON file in which schema is defined in the following manner:

```
import json
 
# Opening JSON file
schema = open('schema_file.json')
 
# returns JSON object as a dictionary
payload = json.load(schema)
```

### insert_schema()
Use `insert_schema(repo_key, payload)` to insert a new schema to an OmixAtlas.

```
omixatlas.insert_schema(repo_key, payload)
```

`repo_key`: (str) repo_id OR repo_name. This is a mandatory field.

`payload`: (dict) The payload is a JSON file which should be as per the structure defined for schema. Only data-admin will have the authentication to update the schema.

`payload` can be loaded from the JSON file in which schema is defined in the following manner:

```
import json
 
# Opening JSON file
schema = open('schema_file.json')
 
# returns JSON object as a dictionary
payload = json.load(schema)
```

### replace_schema()
Replace the existing schema with the new schema of a given OmixAtlas. The old Schema in this case will be completely overridden by the new schema.

```
import json 
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas(token) 
repo_id = "<repo_id>"
with open("<path_of_the_file>/<filename>", "r") as file: 
  schema_data = json.load(file) 

res = omixatlas.replace_schema(repo_id, schema_data) 
```

### validate_schema()

This function is to be used for validating the schema before it’s uploaded to the table of an OmixAtlas

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas (AUTH_TOKEN)
import json

schema = open('schema_file.json') # Opening JSON file

payload = json.load(schema) # returns JSON object as a dictionary

omixatlas.validate_schema(payload)
```

### get_schema()
In order to fetch the schema payload of any existing OmixAtlas, please use get\_schema function of polly-python. Know more in section 4.2.6.2.1 of[ ](https://docs.elucidata.io/OmixAtlas/Polly%20Python.html)[Polly Python - Polly Documentation](https://docs.elucidata.io/OmixAtlas/Polly%20Python.html)

```get_schema(repo_id, schema_level: list (optional), source: str (optional), data_type: str (optional), return_type: (str) (optional)```

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()
schema = omixatlas.get_schema (repo_id,["datasets", "samples", "features"],source="all", data_type="all",return_type="dict")
# to fetch dataset level schema 
schema.datasets
# to fetch sample level schema 
schema.samples
# to fetch feature level schema
schema.features
```

** Preparing the schema in csv file **
A template example of a csv file used to prepare the schema can be found here:
[](https://github.com/ElucidataInc/polly-python/blob/main/Ingest/dataset_schema_case_id_1.csv)[polly-python/dataset_schema_case_id_1.csv at main · ElucidataInc/polly-python](https://github.com/ElucidataInc/polly-python/blob/main/Ingest/dataset_schema_case_id_1.csv)

** Preparing the schema payload to be inserted to an OmixAtlas **

The code to be used to generate the schema payload using the csv file above can be found in this notebook (refer to generating schema payload section) :
[](https://github.com/ElucidataInc/polly-python/blob/main/Ingest/Data_Ingestion_CaseID_1.ipynb)[polly-python/Data_Ingestion_CaseID_1.ipynb at main · ElucidataInc/polly-python](https://github.com/ElucidataInc/polly-python/blob/main/Ingest/Data_Ingestion_CaseID_1.ipynb) 
