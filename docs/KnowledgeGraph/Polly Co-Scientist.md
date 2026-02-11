# Polly Co-Scientist

**Polly Co-Scientist** is an AI-powered research assistant developed to help scientists transform complex biomedical data into actionable scientific insights with ease. It transforms natural language queries into Cypher commands to interact with Polly's Knowledge Graph, allowing users to explore relationships, simulate biological reasoning, and generate hypotheses, without **writing a single line of code.**
By abstracting away query syntax and graph complexity, Co-Scientist allows scientists to focus on scientific reasoning rather than technical execution.  

![KG1](../img/KG/EKG1.png) <center> Polly Co-Scientist Interface</center>

## Working Flow of Polly Co-Scientist

Co-Scientist now collaborates with you to construct queries instead of executing them immediately. Each request follows a structured, five-step workflow designed to ensure accurate intent interpretation and precise Knowledge Graph execution:

- **Intent Parsing**: Automatically extracts key entities, relationships, and filtering criteria from natural language input.

- **Strategic Mapping**: Aligns extracted terms with your Knowledge Graph schema to ensure logical and biologically valid query paths.

- **Output Optimization**: Selects the most appropriate result format—table, sub-graph, or chart—based on the query intent.

- **Cypher Synthesis**: Generates an optimized Cypher query reflecting the approved strategy.

- **User Confirmation (Human-in-the-Loop)**: Presents the complete reasoning and execution plan for review before execution. Users can approve or modify any step, with instant recalibration applied.

- **Execution & Results:** Once confirmed, the query is executed and the results are presented in both textual and visual formats for easy interpretation and exploration.


## Core Enhancements & Optimizations

- **Advanced Relationship Mapping**: Improved recognition of multi-hop relationships, allowing you to discover deep connections between entities.
  
- **Radical Transparency**: You’ll now see a much more detailed breakdown of the internal logic for every step the AI takes.
  
- **Enhanced Consistency**: A set process ensures that the same query yields the same high-quality output every time, regardless of the session.
  
- **Expanded Toolset**: Co-Scientist is now equipped with specialized tools to:
  
    - Extract specific nodes/edges and their counts from the Knowledge Graph (KG).
    
    - Run exploratory queries to understand data distribution.
    
    - Write more efficient, high-performance Cypher queries.


---


## How to Query the Knowledge Graph with Natural Language

Polly Co-Scientist enables you to interact with the biomedical Knowledge Graph using simple, natural language — no technical expertise required. You can run your custom queries. Type your own questions in plain English to retrieve insights tailored to your specific research needs.

### Running Custom Queries on the Knowledge Graph Using Natural Language

#### Step 1: Type Your Query in English
Use the chat interface to enter your research question in natural language.

**Example**:  
`What are the approved drugs for COPD?`

![KG2](../img/KG/EKG2.png) <center> Query</center>

#### Step 2: Review the Interpreted Query Strategy
Co-Scientist analyzes your question and presents a structured interpretation, including:
- Identified research intent

- Suggested entity types and relationship paths

- Related or alternative entities to refine scope

- Recommendations to improve query precision

- If conflicting terms are detected, Co-Scientist prompts you to clarify before proceeding.

You can review, accept, remove, or refine these suggestions before proceeding.

![KG4](../img/KG/EKG4.png) <center> Query Strategy</center>

![KG5](../img/KG/EKG5.png) 


#### Step 3: Inspect the Execution Plan and Generated Cypher Query
After confirmation, Co-Scientist generates and displays the complete execution plan along with the exact Cypher query derived from the approved strategy. This ensures full visibility into how your question will be executed on the Knowledge Graph.

![KG6](../img/KG/EKG6.png) <center> Execution Plan and Cypher Query</center>

#### Step 4: Confirm or Refine (Human-in-the-Loop)
Before execution, Co-Scientist presents its full reasoning and query logic. You can approve the plan as-is or modify any step.
Any change triggers instant recalibration of the execution plan and Cypher query.

![KG7](../img/KG/EKG7.png) 

#### Step 5: View the Results
Once executed, results are presented in the most appropriate format based on the query intent:

**Text Output**
A concise, human-readable summary of the results.

![KG8](../img/KG/EKG8.png) <center> Result</center>


**Visual Output**
An interactive visualization (table, sub-graph, or chart) rendered in the Knowledge Graph viewer for exploration and analysis.

![KG12](../img/KG/EKG12.png) <center> Graph View</center>

![KG9](../img/KG/EKG9.png) <center> Table view</center>





