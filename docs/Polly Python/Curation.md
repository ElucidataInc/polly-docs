::: polly.Curation
    options:
      show_source: false

## Examples

```py
# Install polly python
!sudo pip3 install polly-python --quiet
# Import libraries
from polly.auth import Polly     
from polly.curation  import Curation
import os
import pandas as pd
from json import dumps
import ipywidgets as widgets
```


```py
# Create curation object and authenticate
AUTH_TOKEN=(os.environ['POLLY_REFRESH_TOKEN'])
curate = Curation(AUTH_TOKEN)  
```

### standardize_entity()

```py
# Basic example
curate.standardise_entity("Mouse","species")
```

    {'ontology': 'NCBI',
     'ontology_id': 'txid10090',
     'name': 'Mus musculus',
     'entity_type': 'species',
     'score': None,
     'synonym': None}

```py
# Without 'context'
curate.standardise_entity("AD", "disease")
```

    {'ontology': 'MESH',
     'ontology_id': 'C564330',
     'name': 'Alzheimer Disease, Familial, 3, with Spastic Paraparesis and Apraxia',
     'entity_type': 'disease',
     'score': 202.1661376953125,
     'synonym': 'ad'}

```py
# With context, returns the desired keyword in case of abbreviation
curate.standardise_entity("AD", "disease", 
                context="Patients with atopic dermatitis (AD) where given drug A whereas non AD patients were given drug B")
```

    {'ontology': 'MESH',
     'ontology_id': 'D003876',
     'name': 'Dermatitis, Atopic',
     'entity_type': 'disease',
     'score': 196.61105346679688,
     'synonym': 'atopic dermatitis'}


```py
# Usage of non-matching 'entity_type' returns none values
curate.standardise_entity("Mouse","disease")
```

    {'ontology': 'CUI-less',
     'ontology_id': None,
     'name': None,
     'entity_type': 'disease',
     'score': None,
     'synonym': None}

```py
# Usage of non-supported 'entity_type' returns error -> Here, it is supposed to be "species" and not "specie"
curate.standardise_entity("Mouse","specie")
```

    -----------------------------------------------------------------------

    RequestException                          Traceback (most recent call last)

    Input In [8], in <cell line: 2>()
          1 # Usage of non-supported 'entity_type' returns error -> Here, it is supposed to be "species" and not "specie"
    ----> 2 curate.standardise_entity("Mouse","specie")


    File /usr/local/lib/python3.10/site-packages/polly/curation.py:145, in Curation.standardise_entity(self, mention, entity_type, context, threshold)
        143 if output.get("errors", []):
        144     title, detail = self._handle_errors(output)
    --> 145     raise RequestException(title, detail)
        147 if "term" not in output:
        148     return {
        149         "ontology": "CUI-less",
        150         "ontology_id": None,
        151         "name": None,
        152         "entity_type": entity_type,
        153     }


    RequestException: ('Invalid Payload', {'detail': [{'loc': ['body', 'mention', 'entity_type'], 'msg': "value is not a valid enumeration member; permitted: 'disease', 'drug', 'drug_chebi', 'species', 'tissue', 'cell_type', 'cell_line', 'gene', 'metabolite'", 'type': 'type_error.enum', 'ctx': {'enum_values': ['disease', 'drug', 'drug_chebi', 'species', 'tissue', 'cell_type', 'cell_line', 'gene', 'metabolite']}}]})


### recognise_entity()

```py
# Basic example with two entities
curate.recognise_entity("Gene expression profiling on mice lungs and reveals ACE2 upregulation")
```

    [{'keyword': 'lungs',
      'entity_type': 'tissue',
      'span_begin': 34,
      'span_end': 39,
      'score': 0.9985597729682922},
     {'keyword': 'ACE2',
      'entity_type': 'gene',
      'span_begin': 52,
      'span_end': 55,
      'score': 0.9900580048561096},
     {'keyword': 'mice',
      'entity_type': 'species',
      'span_begin': 29,
      'span_end': 32,
      'score': 0.989605188369751}]

```py
# Multiple entities of the same type
curate.recognise_entity("Batch effects were observed between ductal carcinoma and lobular carcinoma")
```

    [{'keyword': 'ductal carcinoma',
      'entity_type': 'disease',
      'span_begin': 36,
      'span_end': 51,
      'score': 0.9999971389770508},
     {'keyword': 'lobular carcinoma',
      'entity_type': 'disease',
      'span_begin': 57,
      'span_end': 73,
      'score': 0.9999983906745911}]

```py
# Repeating entities
curate.recognise_entity("The study showed ACE2 upregulation and ACE2 downregulation")
```

    [{'keyword': 'ACE2',
      'entity_type': 'gene',
      'span_begin': 17,
      'span_end': 20,
      'score': 0.9962862730026245},
     {'keyword': 'ACE2',
      'entity_type': 'gene',
      'span_begin': 39,
      'span_end': 42,
      'score': 0.990687906742096}]

```py
# No entity in the text
curate.recognise_entity("Significant upregulation was found in 100 samples")
```

    []

### annotate_with_ontology()

```py
# Basic example
curate.annotate_with_ontology("Mouse model shows presence of Adeno carcinoma")
```

    [Tag(name='Mus musculus', ontology_id='NCBI:txid10090', entity_type='species'),
     Tag(name='Adenocarcinoma', ontology_id='MESH:D000230', entity_type='disease')]

```py
# Spelling errors
curate.annotate_with_ontology("Mouse model shows presence of Adino carcinoma")
```

    [Tag(name='Mus musculus', ontology_id='NCBI:txid10090', entity_type='species')]

```py
# incorrect input format -> here, list instead of string
curate.annotate_with_ontology(["Mouse model shows presence", "adeno carcinoma"])
```

    ---------------------------------------------------------------------------

    RequestException                          Traceback (most recent call last)

    Input In [23], in <cell line: 2>()
          1 # incorrect input format -> here, list instead of string
    ----> 2 curate.annotate_with_ontology(["Mouse model shows presence", "adeno carcinoma"])


    File /usr/local/lib/python3.10/site-packages/polly/curation.py:227, in Curation.annotate_with_ontology(self, text)
        209 def annotate_with_ontology(
        210     self,
        211     text: str,
        212 ) -> List[Tag]:
        214     """
        215     Tag a given piece of text. A "tag" is just an ontology term.
        216     Annotates with Polly supported ontologies.
       (...)
        224         tags (set of tuples): set of unique tags
        225     """
    --> 227     entities = self.recognise_entity(text, normalize_output=True)
        228     res = {
        229         self.Tag(
        230             e.get("name", []), e.get("ontology_id", []), e.get("entity_type", [])
       (...)
        233         if e.get("name")
        234     }
        235     return list(res)


    File /usr/local/lib/python3.10/site-packages/polly/curation.py:184, in Curation.recognise_entity(self, text, threshold, normalize_output)
        182 if "errors" in response:
        183     title, detail = self._handle_errors(response)
    --> 184     raise RequestException(title, detail)
        185 try:
        186     entities = response.get("entities", [])


    RequestException: ('Invalid Payload', {'detail': [{'loc': ['body', 'text'], 'msg': 'str type expected', 'type': 'type_error.str'}]})


### find_abbreviations()


```py
# Full form is not mentioned on the text
curate.find_abbreviations("Patient is diagnosed with T1D")
```

    {}

```py
# '-' on the text is not understood
curate.find_abbreviations("Patient is diagnosed with T1D- Type 1 Diabetes")
```

    {}

```py
# Abbreviation is recognized
curate.find_abbreviations("Patient is diagnosed with T1D (Type 1 Diabetes)")
```

    {'T1D': 'Type 1 Diabetes'}

```py
# Abbreviation does not match the full text
curate.find_abbreviations("Patient is diagnosed with T2D (Type 1 Diabetes)")
```

    {}

### assign_control_pert_labels()


```py
sample_metadata = pd.DataFrame({"sample_id": [1, 2, 3, 4], "disease": ["control1", "ctrl2", "healthy", "HCC"],})
sample_metadata
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_id</th>
      <th>disease</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>control1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>ctrl2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>healthy</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>HCC</td>
    </tr>
  </tbody>
</table>
</div>

```py
curate.assign_control_pert_labels(sample_metadata, columns_to_exclude=["sample_id"])
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_id</th>
      <th>disease</th>
      <th>is_control</th>
      <th>control_prob</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>control1</td>
      <td>True</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>ctrl2</td>
      <td>True</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>healthy</td>
      <td>True</td>
      <td>0.96</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>HCC</td>
      <td>False</td>
      <td>0.08</td>
    </tr>
  </tbody>
</table>
</div>


## Tutorial Notebooks

1. [Basic Usage Examples](https://github.com/ElucidataInc/polly-python/blob/main/Curation/Custom%20Curation%20on%20Polly%20Python.ipynb)

2. [Custom Curation with GEO Datasets from Polly](https://github.com/ElucidataInc/polly-python/blob/main/Curation/Using%20the%20Curation%20Library%20on%20Polly-Python.ipynb)
