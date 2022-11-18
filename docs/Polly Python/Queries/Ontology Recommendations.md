Ontology recommendation functionality for disease and tissue are added in Polly-Python. For disease the recommendations are derived from MeSH ontology and for tissue we're using Brenda Tissue Ontology (BTO).

In the existing SQL query itslef, the users would now be able to call a function - 'recommend' on disease and tissue column of metadata to get recommendations.

Usage of 'recommend' function -
```
recommend(field_name, search_term, key - ['match' | 'related'])
```
`field_name`: It can take value: curated_disease, curated_tissue for disease and tissue respectively.

`search_term`: Disease or tissue terms for which recommendations are required.

`key`: Can be "match" which fetches only the terms that have an exact match with the search_term OR "related" which fetches the list of expanded terms synonyms, hypernyms along with match results for the search_term. 

Example:-
```
sql_query = """SELECT dataset_id, curated_disease, curated_tissue FROM geo.datasets WHERE 
        CONTAINS(curated_disease, recommend('curated_disease', 'breast neoplasms', 'related')) AND 
        CONTAINS(curated_tissue, recommend('curated_tissue', 'breast', 'related'))""" 
result = omixatlas.query_metadata(sql_query)
```
For more details and examplpes, please check this [notebook](https://github.com/ElucidataInc/polly-python/blob/main/Discover/ontology_recommendation_disease_tissue.ipynb) 
