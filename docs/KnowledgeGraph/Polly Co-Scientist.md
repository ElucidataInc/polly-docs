# Polly Co-Scientist

**Polly Co-Scientist** is an AI-powered research assistant developed to help scientists transform complex biomedical data into actionable scientific insights with ease. It transforms natural language queries into Cypher commands to interact with Polly's Knowledge Graph, allowing users to explore relationships, simulate biological reasoning, and generate hypotheses, without **writing a single line of code.**

![PollyKG](../img/KG/PollyKG.png) <center> Polly Co-Scientist</center>


## Data Distribution and Sources Used

Polly KG integrates data from multiple open-source and authoritative biomedical databases. Below is a list of sources and their versions:


| DB Name     | Source                                                                                           | Version     |
|-------------|------------------------------------------------------------------------------------------------|-------------|
| MONDO       | [MONDO Ontology](https://mondo.org)                                                            | v2025-02-04 |
| HPO         | [HPO Ontology](https://hpo.jax.org/app/download)                                               | v2025-01-16 |
| NCBI        | [NCBI Download](https://www.ncbi.nlm.nih.gov/home/download/)                                  | release-113 |
| GO          | [Gene Ontology](https://geneontology.org/docs/download-ontology/)                              | v2025-02-06 |
| Reactome    | [Reactome Data](https://reactome.org/download-data)                                            | v91         |
| ChEMBL      | [ChEMBL Downloads](https://chembl.gitbook.io/chembl-interface-documentation/downloads)         | ChEMBL35    |
| OpenTargets | [OpenTargets Platform](https://platform.opentargets.org/downloads) <br> - associationByDatatypeDirect: target-disease associations <br> - diseaseToPhenotype: disease-phenotype mappings <br> - diseases: disease names and IDs <br> - indication: clinical trial indications <br> - interaction: target interactions <br> - targets: target info and tractability | v24.9       |
| UniProt     | [UniProtKB](https://www.uniprot.org/uniprotkb)                                                | 2021_03     |

  <center>Supported Databases <center>


## Key Capabilities

- **Ease of Use**: Designed with a user-friendly interface that simplifies KG exploration. 
- **Natural Language Interface**: Query the Knowledge Graph using simple English, no need to learn Cypher or write Python code.  
- **Automated Query Conversion**: Intelligently translates natural language input into Cypher queries and executes them on the Knowledge Graph.  
- **Dual Output**: View results in both text (summary and Cypher) and graphical (node-edge) formats. 


---


## How to Query the Knowledge Graph with Natural Language

Polly Co-Scientist enables you to interact with the biomedical Knowledge Graph using simple, natural language — no technical expertise required.

You can choose from two options:

- **Use Predefined Templates**: Select from templatized queries available directly within the Polly Co-Scientist interface to quickly explore common biomedical relationships.
- **Run Custom Queries**: Type your own questions in plain English to retrieve insights tailored to your specific research needs.

![PollyKG](../img/KG/predefined.png) <center> Pre-defined questions</center>


### Running Custom Queries on the Knowledge Graph Using Natural Language

#### Step 1: Type Your Query in English
Use the chat interface to enter your research question in natural language.

**Example**:  
`Which tissue express the gene PDE4B, and what are the associated TPM median and cPKG score?`

![Query](../img/KG/customquery.png) <center> Query Example</center>


#### Step 2: Review the Auto-Generated Cypher Query
Polly Co-Scientist automatically translates your natural language input into a Cypher query and JSON output. Both are displayed for reference and executed without requiring manual input.

**Example**:
```cypher
MATCH (t:Tissue)-[:EXPRESSES]->(g:Gene {name: 'PDE4B'})
RETURN t.name AS tissue_name, t.tpm_median AS tpm_median_value
```

![JSON](../img/KG/JSONoutput.png) <center> Cypher Query</center>


### Step 3: View the Results

The output is generated in 2 parts:

- **Text Output:**
A simple readable summary of the query results is displayed.

  Example:
  `The gene PDE4B is expressed in a wide range of tissues including various glands (adrenal, pituitary, prostate), parts of the brain (amygdala, cerebral cortex, hippocampal formation), and major organs such as the heart, liver, and lungs. Additionally, it is found in the thyroid, uterus, bladder, and various regions of the colon and stomach. Unfortunately, the median TPM values for these tissues are not available in the data provided. There is also no cPKG score information available for these tissues.`


- **Graph Output:**
A visual representation of the results is displayed on the Knowledge Graph viewer.

![Output](../img/KG/output.png) <center> Output of the Query</center>

User can click the zoom icon on the right to expand the Knowledge Graph view for easier exploration and full-screen display.

![Output](../img/KG/KGoutput.png) <center> Entire Knowledge Graph Window</center>


