
Curation functions are able to recognise different entities given a text, normalise them based on certain nomenclature such as Polly compatible ontologies. Entities that  are supported are:  "disease", "drug", "species", "tissue", "cell_type", "cell_line", "gene".

[Notebook link](https://github.com/ElucidataInc/PublicAssets/blob/master/internal-user/Polly_notebook_IPC_lib_polly_implementation_dev_test_cases.ipynb)

### annotate_with_ontology ()

Given a text, users can identify and tag entities in a text. Each entity/tag recognised in the text contains the name( word in the text identified), entity_type and the ontology_id.

```
from polly.auth import Polly 
from polly.curation import Curation 
Polly.auth(Token)
obj = Curation()
obj.annotate_with_ontology (text)
#For example
obj.annotate_with_ontology("mus musculus with BRCA gene knocked out") [Tag(name='BRCA1', ontology_id='HGNC: 1100', entity_type= 'gene')]
```

Argument description:-

- text(str): any text or description from which the user wants to identify and tag entities/keywords.

### standardise_entity()
Given a text and the type of entity it is, users can get the Polly compatible ontology for the text such as the MESH ontology. The function also returns a dictionary containing keys and values of the entity type, ontology (such as NCBI, MeSH), ontology ID (such as the MeSH ID), the score (confidence score), and synonyms if any

```
from polly.auth import Polly 
from polly.curation import Curation 
Polly.auth(Token) 
obj = Curation()
obj.standardise_entity(text, entity_type)
#For example 
obj.standardise_entity("Mus musculus","species")

{'ontology': 'NCBI', 'ontology_id': 'txid10090', 'name': 'Mus musculus', 'entity_type': 'species', 'score': None, 'synonym': None}
```

Argument description:-

- text(str): text or word to be standardised.
entity_type(str): the type of entity the given text is such as “species“, “disease“. It can be any one of the supported entity types.
- threshold(int) (Optional Parameter): filter out entities with confidence score less than the threshold. It is recommended to use the default value.
- context(str) (Optional Parameter): text/description to indicate the context in which the text appears. It's used internally for expanding abbreviations.

### recognise_entity()
Users can simply recognise entities (BIOBERT NER model) in a given text without any ontology standardisation (unlike the annotate_with_ontology function above which normalises as well) . A list of spans of identified entities are returned.

```
from polly.auth import Polly
from polly.curation import Curation 
Polly.auth(Token)
obj = Curation()
obj.recognise_entity(text)
#For example
obj.recognise_entity("Adeno carcinoma was observed")
[{"keyword": 'Adeno carcinoma, 'entity_type': 'disease',
'span_begin': 13,
'span_end': 27,
'score': 0.9999943971633911}]
```
