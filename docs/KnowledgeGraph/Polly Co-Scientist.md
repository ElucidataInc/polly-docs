# Polly Co-Scientist

**Polly Co-Scientist** is an AI-powered research assistant developed to help scientists transform complex biomedical data into actionable scientific insights with ease. It transforms natural language queries into Cypher commands to interact with Polly's Knowledge Graph, allowing users to explore relationships, simulate biological reasoning, and generate hypotheses—without **writing a single line of code.**

Grounded in a knowledge graph built from structured datasets and over 37 million PubMed articles, Polly Co-Scientist integrates more than **300 million curated single-cell expression profiles** with high-value sources such as:

- OpenTargets  
- MONDO  
- Reactome  
- ClinVar  
- DepMap  
- ClinicalTrials  
- Therapeutics Data Commons  
- NCBI  
- CHEMBL  

These data sources are enriched with entity-relationship triples to support contextual exploration at tissue and cell-type resolution.


## Key Capabilities

- **Natural Language Interface**: Query the Knowledge Graph using simple English, no need for Cypher or Python expertise.  
- **Automated Query Conversion**: Converts your input into Cypher and runs it on a Neo4j-powered backend.  
- **Dual Output**: View results in both text (summary and Cypher) and graphical (node-edge) formats.  
- **Biological Reasoning**: Explore mechanisms, pathways, and gene/protein relationships at cellular resolution.  
- **Hypothesis Generation**: Identify connections, shortest paths, and co-occurrence between biological entities.  
- **Follow-Up Queries**: Ask follow-ups or refine existing questions using conversational input.  
- **Agentic AI Behavior**: Interacts like a research collaborator—confirming actions, explaining outputs, and adapting to user feedback.


## How to Query the Knowledge Graph with Natural Language

Polly Co-Scientist allows you to query the biomedical Knowledge Graph using plain English. Here's how it works:

### Step 1: Type Your Query in English
Use the chat interface to enter your research question in natural language.

**Example**:  
`Show me the shortest path between JAK1 and NRF1.`


### Step 2: Review the Cypher Query
Polly Co-Scientist converts your query to Cypher without manual intervention and shows it to you for confirmation. Then the query is run.

**Example**:
```cypher
MATCH (p1:Protein {name: "JAK1"}), (p2:Protein {name: "NRF1"})
MATCH path = shortestPath((p1)-[*..10]-(p2))
RETURN path
```

### Step 3: View the Results

The output is generated in 2 parts:

- **Text Output:**
A summary of the relationships between entities is shown in the chat window.

  Example:
  `(:Protein {name: "JAK1"})-[:INTERACTS_WITH]->(:Protein {name: "NodeA"})-[:REGULATES]->(:Protein {name: "NRF1"})`


- **Graph Output:**
A visual representation of the results is displayed on the Knowledge Graph viewer.


### Step 4: Understand the Summary

The system also generates a simple explanation of the output to help you interpret the results.

**Example:** “There are two shortest paths between JAK1 and NRF1 — one via PXDNL and another via CALML5.”


### Step 5: Ask Follow-Up Questions

You can refine or build on your query using natural language.

**Examples:** 
- Only show paths involving CALML5.
- What is the role of CALML5 in this pathway?

>Note: If required, the user can  explicitly ask for the Cypher query as output and ask the bot to edit the query as required in English before asking it to run the query
