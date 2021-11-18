## About Polly Library
Polly Libraries give access to the various capabilities on Polly like querying, filtering and accessing the data on Polly OmixAtlas. It allows access to data in OmixAtlas over any computational platform (like DataBricks, SageMaker, Polly, etc.) of your choice. These functionalities can be accessed through functions in python and [bash](https://docs.elucidata.io/Scaling%20compute/Polly%20CLI%201.html) which can be used over a Terminal.

## About Polly Python
Polly Python library provides convenient access to the above-mentioned functionalities through function in Python language.

## 1 Installation
### 1.1 Install Polly Python using pip

<pre><code>pip install polly-python</code></pre>

## 2 Getting started
### 2.1 Import from libraries

The following libraries need to be imported over the development environment to access the data.

<pre><code>from polly.omixatlas import OmixAtlas
import pandas as pd
from json import dumps</code></pre>

## 3 Authentication
Authentication of the account is required to be able to access the capabilities of the Polly Python library.

### 3.1 Copying authentication token
To get this token, follow the steps below:

1. Go to [Polly](https://polly.elucidata.io)<br>

2. Click the **User Options** icon from the left-most panel<br>

3. Click on **Authentication** on the panel that appears<br>

4. Click on **Copy** to copy the authentication token<br>

### 3.2 Using the token
The following code is required to add the authentication function in the Polly Python library

<pre><code>AUTH_TOKEN = "[authentication_token_copied]"
library_client = OmixAtlas(AUTH_TOKEN)</code></pre>

## 4 Accessing data in OmixAtlas
### 4.1 Calling a function
Use the response from the authentication token to call any function. E.g.
<pre><code>output = library_client.[function()]</code></pre>

The output of the functions is in JSON and/or data frame formats. You can print/download this output.

### 4.2 Functions in Polly Python
#### 4.2.1 Get details of all OmixAtlases
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

#### 4.2.2 Get summary of an OmixAtlas
The following function details a particular OmixAtlas. The **repo_name/repo_id** of this OmixAtlas can be identified by calling the <code>get_all_omixatlas()</code> function.

<pre><code>omixatlas_summary("[repo_id OR repo_name]")</code></pre>
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

#### 4.2.3 Querying the data and the metadata
To access, filter, and search through the metadata schema, the function mentioned below can be used:

<pre><code> query_metadata("[query_written_in_SQL]") </code></pre>
Refer to the Queries section to understand how you could write a query in SQL. The columns returned would depend on the query that was written. The output of the function is a dataframe or a JSON depending on the operations used in the query.

#### 4.2.4 Downloading any dataset
To download any dataset, the following function can be used to get the signed URL of the dataset.

<pre><code> download_data("[repo_name OR repo_id]", "[dataset_id]")</code></pre>

The <code>[repo_name OR repo_id]</code> of this OmixAtlas can be identified by calling the <code>get_all_omixatlas()</code> function. The <code>[dataset_id]</code> can be obtained by querying the metadata at the dataset level using <code>query_metadata("[query written in SQL]")</code>.

The output of this function is a *signed URL*. The data can be downloaded by clicking on this URL.

> **_NOTE:_** This signed URL expires after 60 minutes from when it is generated.

<br>The output data is in .gct/h5ad format. This data can be parsed into a data frame for better accessibility using the following code:

##### 4.2.4.1 Downloading .gct as a data frame
<pre><code>url = library_client.download_data("[repo_id OR repo_name]", "[dataset_id]").get('data')
file_name = "[dataset_id].gct"
os.system(f"wget -O '{file_name}' '{url['attributes']['download_url']}'")
from cmapPy.pandasGEXpress.parse_gct import parse
data = parse(file_name)</code></pre>

##### 4.2.4.2 Downloading h5ad as a data frame
<pre><code>url = library_client.download_data("[repo_id OR repo_name]", "[dataset_id]").get('data')
file_name = "[dataset_id].h5ad"
os.system(f"wget -O '{file_name}' '{url['attributes']['download_url']}'")
import scanpy as sc
data = sc.read_h5ad(file_name)</code></pre>

#### 4.2.5 Working with workspaces
Polly python provides two functions to create new workspaces and listing existing workspaces
Following library need to be imported to work with workspaces

```
from polly import new_workspaces

library_client = new_workspaces.Workspaces(AUTH_TOKEN)
```
`AUTH_TOKEN` here is the authentication token described in section 3.1 and 3.2

##### 4.2.5.1 Creating new workspace
use `create_workspace` to create a new workspace with desired name
```
library_client.create_workspace("name_of_workspace")
```
##### 4.2.5.2 Fetching existing workspaces
use `fetch_my_workspaces()` to fetch existing workspaces
```
library_client.fetch_my_workspaces()
```  


### 4.3 Data Schema
The data available within OmixAtlas is curated within 5 indexes/tables on the basis of the information it contains. These five indexes are:

**Dataset level metadata (index: files)**: Contains curated fields like drug, disease, tissue organism, etc for each dataset.

**Sample level metadata (index: gct_metadata)**: As the name suggests, this contains sample level metadata information for all samples except single-cell samples. It contains curated fields like cell lines, experimental design, etc for each sample except single-cell samples.

**Sample level metadata for single-cell (index: h5ad_metadata)**: This table only contains metadata information for single-cell samples. It contains curated fields like cell line, gene counts, UMI counts for each sample.

**Feature level metadata (index: gct_data)**: This table contains feature level metadata information for all data except single-cell.  It contains the gene/molecule symbol along with the feature intensity for each sample.

**Feature level metadata for single-cell (index: h5ad_data)**: This table contains feature level metadata information for all single-cell data.  It contains the gene/molecule symbol studied in an experiment along with the cell type and expression value.

To find relevant information that can be used for querying, refer the curated data schema [here](https://docs.elucidata.io/OmixAtlas/Data%20Schema.html).

### 4.4 Writing a query

> **NOTE:** Two API versions, v1 and v2, are currently available in polly-python which refer to infrastructure used for query. This could be specified by setting `query_api_version` to 'v1' or 'v2' while passing the query to a function. For example, `query_api_version='v2'` makes the query against version 2. By default, queries are made against v2. However, not all features are compatible with V2. This document specifies the version compatibility while describing a function or feature.

The complete syntax for searching and aggregating data is as follows:

<pre><code>SELECT [DISTINCT] (* | expression) [[AS] alias] [, ...]
[ WITH with_query [, ...] ]
SELECT [ ALL | DISTINCT ] select_expression [, ...]
[ FROM from_item [, ...] ]
[ WHERE condition ]
[ GROUP BY [ ALL | DISTINCT ] grouping_element [, ...] ]
[ HAVING condition ]
[ { UNION | INTERSECT | EXCEPT } [ ALL | DISTINCT ] select ]
[ ORDER BY expression [ ASC | DESC ] [ NULLS FIRST | NULLS LAST] [, ...] ]
[ OFFSET count [ ROW | ROWS ] ]
[ LIMIT [ count | ALL ] ]
</code></pre>


#### 4.4.1 Querying the dataset level metadata:
<pre><code> query = "SELECT [column_name] FROM [files] WHERE [column_name]='[value]'"
 query = "SELECT [column_name] FROM [repo_name].datasets WHERE [column_name]='[value]'"
 query = "SELECT [column_name] FROM [repo_id].datasets WHERE [column_name]='[value]'"
</code></pre>

#### 4.4.2 Querying the sample level metadata:
##### For all samples except Single Cell
<pre><code>query = "SELECT [column_name] FROM [gct_metadata] WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_name].samples WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_id].samples WHERE [column_name]='[value]'"
</code></pre>

##### For samples in Single Cell
<pre><code>query = "SELECT [column_name] FROM [h5ad_metadata] WHERE [column_name]='[value]'"</code></pre>

#### 4.4.3 Querying the feature level metadata:
##### For all features except Single Cell
<pre><code>query = "SELECT [column_name] FROM [gct_data] WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_name].features WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_id].features WHERE [column_name]='[value]'"
</code></pre>

##### For features in Single Cell
<pre><code>query = "SELECT [column_name] FROM [h5ad_data] WHERE [column_name]='[value]'"</code></pre>

#### 4.4.4 Querying available tables in an altas

The tables available for querying in an atlas could be shown by using the following command

```
query = "SHOW TABLES IN repo_name "
```
For example, all the tables in `liveromix_atlas` could be known by constructing following query and passing it to 'query_metadata'

```
query = "SHOW TABLES IN liveromix_atlas "
library_client.query_metadata(
    query,
    query_api_version="v2"  
)

```
#### 4.4.5 Description of a table

To obtain the names of columns in a table, use following query

```
query = "DESCRIBE repo_name.table_name"
```
For example, to know the dataset columns in `liveromix_atlas`, following code can be used

```
query = " DESCRIBE liveromix_atlas.datasets "
library_client.query_metadata(
    query,
    query_api_version="v2"    
)

```  
#### 4.4.6 Querying expression matrices of datasets

The expression matrices of the datasets can be queried. To do this, the value of `kw_data_table_name` column for a dataset can be used in conjunction with the `expr_tables` prefix.

For example, to search dataset “GSE54234_GPL13912” in liveromix_atlas, the query will look like

```
SELECT * FROM expr_tables.liveromix_atlas__gse54234_gpl13912 WHERE ...

```  
For getting the columns to search in the above example, a query can be constructed to describe `expr_tables.liveromix_atlas__gse54234_gpl13912` table as explained in section 4.4.5

### 4.5 Writing conditions with operators
The following operators can be used to define the conditions in the above mentioned queries:

Operators  | Functions performed
------------- | -------------
<code>=</code> |  **Equal to** operator which can be used to find matching strings with values in the columns
<code><></code> | **Not equal to** operator which can be used to find non-matching strings with values in the columns
<code>></code> | **Greater than** operator which can be used **ONLY** for integer based columns
<code><</code> | **Less than** operator which can be used **ONLY** for integer based columns
<code>>=</code> | **Greater than or equal to** operator which can be used **ONLY** for integer based columns
<code><=</code> | **Less than or equal to** operator which can be used **ONLY** for integer based columns
<code>IS NULL</code> | Check if the field value is <code>NULL</code>.
<code>IS NOT NULL</code> | Check if the field value is <code>NOT NULL</code>.
<code>AND</code> | All values across the parameters searched for have to be present in a dataset for it to be returned as a match when the AND operator is used. <br><br>e.g. “organism = ‘Homo sapiens' AND disease = 'Carcinoma, Hepatocellular’” would only return datasets that belong to homo sapiens and have the disease as hepatocellular carcinoma.
<code>OR</code> | Atleast any one value across the parameters searched for have to be present in a dataset for it to be returned as a match when the OR operator is used. <br><br>e.g. <code>organism = 'Homo sapiens' OR disease = 'Carcinoma, Hepatocellular'</code> would return datasets that belong to homo sapiens or have the disease as hepatocellular carcinoma or match both criteria.
<code>MATCH QUERY(<column_name>,'value')</code> | It works like a fuzzy search. If you add a string for a parameter with this operator, it would return all possible results matching each word in the string. The search output is returned with a “Score” using which the output is sorted. <br><br>e.g. <code>MATCH_QUERY(description,'Transcriptomics profiling')</code> would return all datasets having <code>transcriptomics profiling</code> , <code>Transcriptomics</code> and <code>profiling</code> as possible terms within their description. Each dataset would be scored on the basis of matching of the searched string with the information present within the dataset. It works with API version V1 only hence `query_api_version="v1"` must be passed along with the query.
<code>MATCH PHRASE(<column_name>,'value')</code> | This can be used for exact phrase match with the information being searched for. <br><br>e.g. <code>MATCH_PHRASE(description,'Transcriptomics profiling')</code> would only return the datasets that have <code>Transcriptomics profiling</code> within their description. It works with API version V1 only hence `query_api_version="v1"` must be passed along with the query.
<code>MULTI MATCH('query'='value', 'column_name'='value)</code> | This can be used to search for text in multiple fields, use <code>MULTI MATCH('query'='value', 'column_name'='value)</code>. <br><br>e.g. <code>MULTI MATCH('query'='Stem Cells', 'fields'='tissue','description')</code> would return datasets that have <code>"Stem Cell"</code> in either <code>tissue</code> OR <code>description</code> fields.
<code>GROUP BY</code> | The <code>GROUP BY</code> operator groups rows that have the same values into summary rows. The GROUP BY statement is often used with aggregate functions (COUNT, MAX, MIN, SUM, AVG) to group the result-set by one or more columns.
<code>HAVING</code> | Use the HAVING clause to aggregate inside each bucket based on aggregation functions (COUNT, AVG, SUM, MIN, and MAX). The HAVING clause filters results from the GROUP BY clause
<code>COUNT(*)</code> | This counts each row present in a table/index being queried.<br><br> **NOTE: The output of this query would return a JSON stating the total number of rows in the table**
<code>LIMIT</code> | **NOTE: The response of any query returns 200 entries by default**. <br>You can extend this by defining the LIMIT of the results you want to query to be able to return.
<code>ORDER BY</code> | Can only be used to sort the search results using integer based parameters in the schema. Sorting on the basis of dataset_id, number of samples, <code>_score</code> of the data is available at the dataset-level metadata. <code>ASC</code> or <code>DESC</code> can be used to define whether you want to order the rows in ascending or descending order respectively

## 5 Experimental Expanding Search Feature
A disease is described by many different names in the literature. If we search it with a keyword, datasets annotated with related keywords for the same disease are missed in the output. For example, if 'nash' is searched for disease, 'non-alcoholic fatty liver disease', 'nonalcoholic steatohepatitis', 'nash-non-alcoholic steatohepatitis', 'non-alcoholic steatohepatitis' are missed although they are valid hits for disease keyword 'nash'.

To address this, we have included an experimental feature which allows expansion of disease query to include additional similar keywords based on ontology mapping. The query written in the Polly Python Library is expanded on the basis of the ontology tree mapping available in the following resources

1. [The disease ontology](https://disease-ontology.org/)

2. [Experimental factor ontology](https://www.ebi.ac.uk/efo/)

3. [Mondo Disease Ontology](http://obofoundry.org/ontology/mondo.html)

4. [NCI Thesaurus](https://github.com/NCI-Thesaurus/thesaurus-obo-edition)

5. [Human phenotype ontology](https://hpo.jax.org/)

6. [Ontology for MIRNA Target](https://github.com/OmniSearch/omit)  

Expanding search feature could be use in the following way

> Note: Experimental search feature currently only works with API v1. For this reason, `query_api_version="v1"` must be included while using this feature.    

1. Setting "expand" to `True`: It allows expansion of the query to include disease terms children from the ontology tree. It means no false positive datasets are included in the output while including more true positives.  

```
result = library_client.query_metadata(
    query=user_query,
    experimental_features = {"expand":True, "related_terms":False},
    query_api_version="v1"
  )
```
2. Setting "related_terms" to `True`: It allows expansion of query by including immediate parent of the disease term in the ontology tree. It causes maximum expansion but a few false positive datasets may be included.

```
result = library_client.query_metadata(
    query=user_query,
    experimental_features = {"expand":True, "related_terms":True},
    query_api_version="v1")
```

## 6 Example Use Cases

### 6.1 Querying datasets in Liver OmixAtlas
1. To identify datasets belonging to the tissue Liver, disease Liver cirrhosis and organism Homo sapiens

    ```
        query = """SELECT * FROM liveromix_atlas.datasets
                            WHERE disease = 'liver cirrhosis'
                            AND tissue = "liver"
                            AND organism="Homo sapiens"
                            """
    ```

2. To identify all datasets belonging to a Hepatocellular Carcinoma disease in Human and Mouse

    ```
        query = """SELECT * FROM liveromix_atlas.datasets
                            WHERE disease = 'Carcinoma, Hepatocellular'
                            AND ( organism="Homo sapiens" or organism="Mus musculus")
                          """
    ```

3. To identify type of data, tissue, and data description of all the datasets belonging to a Hepatocellular Carcinoma disease in Human and Mouse

    ```
        query = """SELECT dataset_id,description,tissue,kw_data_type
                            FROM liveromix_atlas.datasets
                            WHERE disease = 'Carcinoma, Hepatocellular'
                            AND ( organism="Homo sapiens" or organism="Mus musculus")
                          """
    ```

4. Identify all transcriptome datasets in Hepatocellular Carcinoma disease in Human and Mouse

    ```
        query = """SELECT * FROM liveromix_atlas.datasets
                        WHERE disease = 'Carcinoma, Hepatocellular'
                        AND ( organism="Homo sapiens" or organism="Mus musculus")
                        AND ( kw_data_type = "Transcriptomics")
                      """
    ```

5. Identify all transcriptome datasets from GEO in Hepatocellular Carcinoma disease in Human and Mouse

    ```
        query = """SELECT * FROM liveromix_atlas.datasets
                        WHERE disease = 'Carcinoma, Hepatocellular'
                        AND ( organism="Homo sapiens" or organism="Mus musculus")
                        AND ( kw_data_type = "Transcriptomics")
                        AND ( dataset_source = "GEO")
                      """
    ```

6. Select all methylation, miRNA, and transcriptomics datasets for Hepatocellular Carcinoma disease in Human and Mouse

    ```
        query = """SELECT * FROM liveromix_atlas.datasets
                        WHERE dataset_source = 'TCGA'
                        AND (kw_data_type = 'Methylation'
                            OR kw_data_type = 'MiRNA'
                            OR kw_data_type = 'Transcriptomics'
                            )
                        AND ( organism="Homo sapiens" or organism="Mus musculus")
                      """
    ```
In all the examples above `liveromix_atlas_files` could be used in place of `liveromix_atlas.datasets`

### 6.2 Querying samples in Liver OmixAtlas
1. Select a sample with a GEO ID

    ```
        query = """SELECT *
                FROM liveromix_atlas.samples
                WHERE kw_column  = 'GSM798352'
                """
    ```

2. Identify all samples in which "CYP1B1" gene has been knocked out

    ```
        query = """SELECT *
                FROM liveromix_atlas.samples
                WHERE kw_curated_genetic_mod_type = 'knockout'
                AND kw_curated_modified_gene = "CYP1B1"
                """
    ```

3. Identify all samples where gene "YAP1" has been knocked in, in hepatocype cell line

    ```
        query = """SELECT kw_curated_cell_line,kw_curated_drug,kw_curated_genetic_mod_type,kw_curated_modified_gene,kw_curated_cell_type
                    FROM liveromix_atlas.samples
                    WHERE kw_curated_genetic_mod_type = 'knockin'
                    AND kw_curated_cell_type = "hepatocyte"
                    AND kw_curated_modified_gene = "YAP1"
                """
    ```

4. Identify all samples of 'Carcinoma, Hepatocellular' disease that have been treated with some drug

    ```
        query = """SELECT kw_curated_cell_line,kw_curated_drug,kw_curated_tissue,kw_curated_disease
                FROM liveromix_atlas.samples
                WHERE kw_curated_disease = 'Carcinoma, Hepatocellular'
                AND kw_curated_drug <> "none"
                """
    ```
In all the examples above `liveromix_atlas_gct_metadata` could be used in place of `liveromix_atlas.samples`

### 6.3 Querying features in Liver OmixAtlas
1. Identify all samples in which the gene “YAP1” is upregulated

    ```
        query = """SELECT *
                    FROM liveromix_atlas.features
                    WHERE kw_index = 'YAP1'
                    AND kw_column.kw_expression > 0
                """
    ```
In this example `liveromix_atlas_gct_data` could be used in place of `liveromix_atlas.features`
