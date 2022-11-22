### get_all_omixatlas()

Get details of all the OmixAtlases accessible by you.

<pre><code>omixatlas.get_all_omixatlas() </code></pre>

The output of this function would be JSON containing
<pre><code>{'data':
[
{
  "repo_name": "geo",
  "repo_id": "9",
  "indexes": {
    "gct_metadata": "geo_gct_metadata",
    "h5ad_metadata": "geo_h5ad_metadata",
    "csv": "geo_csv",
    "files": "geo_files",
    "ipynb": "geo_ipynb",
    "gct_data": "geo_gct_data",
    "h5ad_data": "geo_h5ad_data"
  },
  "v2_indexes": {
    "csv": "geo_csv",
    "files": "geo_files",
    "gct_col_metadata": "geo_gct_col_metadata",
    "gct_row_metadata": "geo_gct_row_metadata",
    "h5ad_col_metadata": "geo_h5ad_col_metadata",
    "h5ad_row_metadata": "geo_h5ad_row_metadata",
    "ipynb": "geo_ipynb",
    "json": "geo_json"
  },
  "sources": [
    {
      "geo": 98902
    }
  ],
  "datatypes": [
    {
      "transcriptomics": 98863
    },
    {
      "raw counts transcriptomics": 38
    }
  ],
  "dataset_count": 98902,
  "disease_count": 3760,
  "tissue_count": 1078,
  "organism_count": 317,
  "cell_line_count": 4492,
  "cell_type_count": 904,
  "drug_count": 1394,
  "data_type_count": 2,
  "data_source_count": 1,
  "sample_count": 2286352,
  "normal_sample_count": 1575120
}
    {...},
    {...}
]
}</code></pre>

### omixatlas_summary()

Get summary of a particular OmixAtlas. The **repo_name/repo_id** of this OmixAtlas can be identified by calling the <code>get_all_omixatlas()</code> function.

<pre><code>omixatlas.omixatlas_summary("[repo_id OR repo_name]")</code></pre>
The output of this function would be JSON containing

<pre><code>
{
  "data": {
    "repo_name": "tcga",
    "repo_id": "15",
    "indexes": {
      "gct_metadata": "tcga_gct_metadata",
      "h5ad_metadata": "tcga_h5ad_metadata",
      "csv": "tcga_csv",
      "files": "tcga_files",
      "ipynb": "tcga_ipynb",
      "gct_data": "tcga_gct_data",
      "h5ad_data": "tcga_h5ad_data"
    },
    "v2_indexes": {
      "csv": "tcga_csv",
      "files": "tcga_files",
      "gct_col_metadata": "tcga_gct_col_metadata",
      "gct_row_metadata": "tcga_gct_row_metadata",
      "h5ad_col_metadata": "tcga_h5ad_col_metadata",
      "h5ad_row_metadata": "tcga_h5ad_row_metadata",
      "ipynb": "tcga_ipynb",
      "json": "tcga_json"
    },
    "sources": [
      {
        "tcga": 55062
      }
    ],
    "datatypes": [
      {
        "transcriptomics": 11093
      },
      {
        "mirna": 10468
      },
      {
        "copy number variation": 10411
      },
      {
        "mutation": 8810
      },
      {
        "methylation": 8027
      },
      {
        "proteomics": 6253
      }
    ],
    "dataset_count": 55062,
    "disease_count": 34,
    "tissue_count": 58,
    "organism_count": 1,
    "cell_line_count": 0,
    "cell_type_count": 0,
    "drug_count": 812,
    "data_type_count": 6,
    "data_source_count": 1,
    "sample_count": 47143,
    "normal_sample_count": 0
  }
}
</code></pre>

### create()

Data-admins can create an Omixatlas using polly-python. The function create takes in four parameters as described below.
<pre><code>from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()
new_repo = omixatlas.create("[display_name]", "[description]", 
                          repo_name ="[repo_name]" (optional), 
                          image_url = "[image_url]" (optional))</code></pre>

Constraints on the parameters:-

`display_name` (str) Alphanumeric characters are allowed and the length constraint is between 1 to 30 characters.

`description` (str) All characters are allowed and the length constraint is between 1 to 100 characters.

`image_url` (str) Users can also enter the path of image_url that they want to showcase on the newly created Omixatlas tile. If the image_url is not provided then the system puts up a default image on the tile of the newly created Omixatlas. [example_image](https://elucidatainc.github.io/PublicAssets/discover-fe-assets/omixatlas_hex.svg)

`repo_name` (str) Lowercase alphanumeric characters (separated by _) is allowed and between 1 to 30 characters.

### update()

Data-admin can update the following metadata of an Omixatlas:-

a. Metadata for an existing Omixatlas. The attributes that can be updated are display_name, description, and image_url.

b. Adding components(i.e. Apps, Python Notebooks) to an existing Omixatlas.

<pre><code>from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas()
omixatlas.update(repo_key, display_name = <Updated display_name string>, description = <Updated description string>,)
                image_url = <image_url string>, components = [component_1])
                
component_1 = {"data_type":["string"], "component_id":integer}
// example:- component_1 = {"data_type":["Transcriptomics"], "component_id":78}</code></pre>

`repo_key` (str or int) Repository ID or repository name of the OmixAtlas 

`display_name` (str) The name to be displayed on the front-end

`description` (str) Description for the OmixAtlas

`image_url` (str) The logo of the OmixAtlas as shown on the front-end. The link provided here should be hosted on a public domain.

`components` Application or notebooks which should be linked with given data types of an OmixAtlas

Constraints on components: Components will be a dictionary that will have two mandatory key-value pairs that are data_type and component_id.
