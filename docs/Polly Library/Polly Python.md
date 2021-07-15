## About Polly Library
Polly Libraries give access to the various capabilities on Polly like querying, filtering and accessing the data on Polly OmixAtlas. It allows access to data in OmixAtlas over any computational platform (like DataBricks, SageMaker, Poly, etc.) of your choice. These functionalities can be accessed through functions in python and [bash](https://docs.elucidata.io/Scaling%20compute/Polly%20CLI%201.html)

## About Polly Python 
Polly Python library provides convenient access to the above-mentioned functionalities through function in Python language.

## Installation
### Install Polly Python using pip

<pre><code>pip install [version] </code></pre>

## Getting started
### Import from libraries

The following libraries need to be imported over the development environment to access the data.

<pre><code>from polly.omixatlas import OmixAtlas
import pandas as pd
from json import dumps</code></pre>

## Authentication
Authentication of the account is required to be able to access the capabilities of the Polly Python library.

### Copying the token for authentication
1. Go to [Polly](https://polly.elucidata.io)

2. Click the *User Options* icon from the left-most panel

3. Click on *Authentication* on the panel that appears

4. Click on *Copy* to copy the authentication token

### Using the token
The following code is required to add the authentication function in the Polly Python library

<pre><code>AUTH_TOKEN = "[authentication_token_copied]"
library_client = OmixAtlas(AUTH_TOKEN)</code></pre>

## OmixAtlas
### Calling a function
Use the response from the authentication token to call any function. E.g.
<pre><code>output = library_client.[function()]</code></pre>

The output of the functions is in JSON and/or data frame formats. You can print/download this output.

### Functions in Polly Python
#### 1. Get details of all OmixAtlases
The following function details all the OmixAtlases accessible by you.

<pre><code> get_all_omixatlas() </code></pre>

The output of this function would be JSON containing
<pre><code>{'data': 
[
  {'repo_name': 'name', 
    'repo_id': 'id', 
    'indexes': 
    { 
      'gct_metadata': 'abc', 
      'h5ad_metadata': 'abc', 
      'csv': 'abc', 
      'files': 'abc', 
      'json': 'abc', 
      'ipynb': 'abc', 
      'gct_data': 'abc', 
      'h5ad_data': 'abc'
    }, 
    'dataset_count': 123, 
    'disease_count': 123, 
    'diseases': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'organism_count': 123, 
    'organisms': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'sources': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'datatypes': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'sample_count': 123
    }, 
    {...},
    {...}
  ]
}</code></pre>

#### 2. Get the summary of any OmixAtlas
The following function details a particular OmixAtlas. The key/repo id of this OmixAtlas can be identified by calling the get_all_omixatlas() function.

<pre><code>omixatlas_summary(”[repo_id OR repo_name]”)</code></pre>
The output of this function would be JSON containing

<pre><code>{'data': 
  {
    'repo_name': 'name', 
    'repo_id': 'id', 
    'indexes': 
    { 
      'gct_metadata': 'abc', 
      'h5ad_metadata': 'abc', 
      'csv': 'abc', 
      'files': 'abc', 
      'json': 'abc', 
      'ipynb': 'abc', 
      'gct_data': 'abc', 
      'h5ad_data': 'abc'
    }, 
    'dataset_count': 123, 
    'disease_count': 123, 
    'diseases': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'organism_count': 123, 
    'organisms': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'sources': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'datatypes': ['abc', 'bcd', 'cde', 'def', 'efg', 'fgh', 'ghi', 'hij', 'ijk'], 
    'sample_count': 123
  }
}</code></pre>

#### 3. Querying the data and the metadata
To access, filter, and search through the metadata schema, the function mentioned below can be used:

<pre><code> query_metadata(“[query_written_in_SQL]”) </code></pre>
Refer to the Queries section to understand how you could write a query in SQL. The columns returned would depend on the query that was written. The output of the function is a dataframe or a JSON depending on the operations used in the query. 

#### 4. Downloading any dataset
To download any dataset, the following function can be used to get the signed URL of the dataset. The key/repo id of this OmixAtlas can be identified by calling the get_all_omixatlas() function. The _id can be obtained from querying the metadata at the dataset level using query_metadata(“<query written in SQL>”) to get the _id column.

<pre><code> download_data(”[repo_name OR repo_id]”, “[value_of _id_column_for_a_dataset]”)</code></pre>

The output of this function is a *signed URL*. The data from this URL can be downloaded using requests library/other libraries in python.

E.g. Download the data from the URL as a .csv using *requests* library

<pre><code>import requests
url = "example.com"

r = requests.get(url)  
with open("[name_of_file].csv",'wb') as f:
    f.write(r.content)</code></pre>
  
### Queries
#### The syntax for querying the dataset level metadata:
<pre><code> query = “SELECT [column_name] FROM liver_atlas_files WHERE [column_name]='[value]’” </code></pre>

####The syntax for querying the sample level metadata:
##### For all samples except Single Cell
<pre><code>query = “SELECT [column_name] FROM liver_atlas_gct_metadata WHERE [column_name]='[value]’”</code></pre>

##### For samples in Single Cell
<pre><code>query = “SELECT [column_name] FROM liver_atlas_h5ad_metadata WHERE [column_name]='[value]’”</code></pre>

#### The syntax for querying the feature level metadata:
##### For all features except Single Cell
<pre><code>query = “SELECT [column_name] FROM liver_atlas_gct_data WHERE [column_name]='[value]’”</code></pre>

##### For features in Single Cell
<pre><code>query = “SELECT [column_name] FROM liver_atlas_h5ad_data WHERE [column_name]='[value]’”</code></pre>

### Operators

Operators  | Functions performed | Ouput
------------- | ------------- | ------------
= |  *Equal to* operator which can be used to find matching strings with values in the columns | DataFrame
!= | *Not equal to* operator which can be used to non-matching strings with values in the columns | DataFrame
