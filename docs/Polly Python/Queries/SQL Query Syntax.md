#### Complete Syntax:
The complete syntax for searching and aggregating data is as follows:

<pre><code>[ WITH with_query [, ...] ]
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

#### Querying the dataset level metadata:
<pre><code> query = "SELECT [column_name] FROM [files] WHERE [column_name]='[value]'"
 query = "SELECT [column_name] FROM [repo_name].datasets WHERE [column_name]='[value]'"
 query = "SELECT [column_name] FROM [repo_id].datasets WHERE [column_name]='[value]'"
</code></pre>

#### Querying the sample level metadata:
##### For all samples except Single Cell
<pre><code>query = "SELECT [column_name] FROM [gct_metadata] WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_name].samples WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_id].samples WHERE [column_name]='[value]'"
</code></pre>

##### For samples in Single Cell
<pre><code>query = "SELECT [column_name] FROM [h5ad_metadata] WHERE [column_name]='[value]'"</code></pre>

#### Querying the feature level metadata:
##### For all features except Single Cell
<pre><code>query = "SELECT [column_name] FROM [gct_data] WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_name].features WHERE [column_name]='[value]'"
query = "SELECT [column_name] FROM [repo_id].features WHERE [column_name]='[value]'"
</code></pre>

##### For features in Single Cell
<pre><code>query = "SELECT [column_name] FROM [h5ad_data] WHERE [column_name]='[value]'"</code></pre>

#### Query specific to source and datatype in an OmixAtlas

If the schema of the repository has multiple source and data types then querying for the specific source and/or datatype is enabled in this release.

```
from polly import Omixatlas
omixatlas = Omixatlas()
query = """SELECT * FROM repo_name.source_name_in_schema.datatype_name_in_schema.datasets WHERE CONTAINS(curated_disease, 'Multiple Myeloma')
results-omixatlas.query_metadata(query)

results
```

For clearer understanding of this feature, let’s look at a case:-

In the schema shown below for repo_id: 1659450268526, we can see there are 3 different sources. Out of these sources, the source lincs has 2 datatypes - mutation and transcriptomics

![](../img/OmixAtlas-Images/29_dataset_level_schema.png)

<p align=center>dataset level schema of repo_id 1659450268526
</p>

Query examples which will be supported in the above context:-

1 query = """SELECT * FROM 1659450268526.geo.datasets WHERE CONTAINS(curated_disease, 'Multiple Myeloma') """

2 query = """SELECT * FROM 1659450268526.lincs.datasets WHERE CONTAINS(curated_disease, 'Multiple Myeloma') """

3 query = """SELECT * FROM 1659450268526.lincs.transcriptomics.datasets WHERE CONTAINS(curated_disease, 'Multiple Myeloma') """


The response to these queries will have a dataframe with columns which are specific to the source and datatype as per the schema shown above.  

A detailed notebook and technical note describing these queries will be updated in the documentation.

#### Writing conditions with operators
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
<code>MATCH QUERY(<column_name>,'value')</code> | It works like a fuzzy search. If you add a string for a parameter with this operator, it would return all possible results matching each word in the string. The search output is returned with a “Score” using which the output is sorted. <br><br>e.g. <code>MATCH_QUERY(description,'Transcriptomics profiling')</code> would return all datasets having <code>transcriptomics profiling</code> , <code>Transcriptomics</code> and <code>profiling</code> as possible terms within their description. Each dataset would be scored on the basis of matching of the searched string with the information present within the dataset.
<code>MATCH PHRASE(<column_name>,'value')</code> | This can be used for exact phrase match with the information being searched for. <br><br>e.g. <code>MATCH_PHRASE(description,'Transcriptomics profiling')</code> would only return the datasets that have <code>Transcriptomics profiling</code> within their description.
<code>MULTI MATCH('query'='value', 'column_name'='value)</code> | This can be used to search for text in multiple fields, use <code>MULTI MATCH('query'='value', 'column_name'='value)</code>. <br><br>e.g. <code>MULTI MATCH('query'='Stem Cells', 'fields'='tissue','description')</code> would return datasets that have <code>"Stem Cell"</code> in either <code>tissue</code> OR <code>description</code> fields.
<code>GROUP BY</code> | The <code>GROUP BY</code> operator groups rows that have the same values into summary rows. The GROUP BY statement is often used with aggregate functions (COUNT, MAX, MIN, SUM, AVG) to group the result-set by one or more columns.
<code>HAVING</code> | Use the HAVING clause to aggregate inside each bucket based on aggregation functions (COUNT, AVG, SUM, MIN, and MAX). The HAVING clause filters results from the GROUP BY clause
<code>COUNT(*)</code> | This counts each row present in a table/index being queried.<br><br> **NOTE: The output of this query would return a JSON stating the total number of rows in the table**
<code>LIMIT</code> | **NOTE: The response of any query returns 200 entries by default**. <br>You can extend this by defining the LIMIT of the results you want to query to be able to return.
<code>ORDER BY</code> | Can only be used to sort the search results using integer based parameters in the schema. Sorting on the basis of dataset_id, number of samples, <code>_score</code> of the data is available at the dataset-level metadata. <code>ASC</code> or <code>DESC</code> can be used to define whether you want to order the rows in ascending or descending order respectively



