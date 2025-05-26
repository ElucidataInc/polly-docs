# Terminology in the Polly KG Explorer Application

This section provides definitions of key terms used in the Polly KG Explorer Application. Understanding these terms will help users navigate and utilize the tool more effectively.

## General Terms

- **Knowledge Graph (KG)**: A structured representation of relationships between biological entities such as genes, proteins, diseases, and drugs. It enables data integration and discovery of new insights.
- **Node**: A fundamental unit in the knowledge graph representing a biological entity (e.g., a gene, protein, disease, or drug).
- **Edge**: A connection between two nodes that represents a relationship between them (e.g., a gene encodes a protein).
- **Subgraph**: A smaller, extracted portion of the knowledge graph that focuses on selected nodes and their relationships.
- **Hops**: The number of steps (edges) between two nodes in the graph. A 1-hop relationship means nodes are directly connected, while a 2-hop relationship means they are connected through an intermediate node.
- **Searched Node**: The primary node that the user searches for in the application. This node serves as the central point for exploring relationships within the knowledge graph.
- **Immediate Neighbors**: Nodes that share a direct edge with the Searched Node. These nodes represent entities that have a direct relationship with the searched entity (e.g., a drug directly linked to a target protein).

## Application-Specific Terms

- **Node Selection Panel**: The section of the interface where users can select node categories and specific values for visualization.
- **Subgraph Visualization Area**: The main workspace where the selected nodes and their relationships are displayed as an interactive graph.
- **Graph Legend**: A reference guide within the application that explains the color coding for different node categories.
- **2D/3D Toggle**: A feature allowing users to switch between two-dimensional and three-dimensional visual representations of the knowledge graph.
- **Generate Button**: The action button that creates or updates the subgraph based on the selected nodes.

## Data & Case Sensitivity Terms

**Case Sensitivity in Node Search**: The application requires specific casing for different node types:

- **Drugs**: All uppercase (e.g., `ASPIRIN`)
- **Genes**: All uppercase (e.g., `ALPK1`)
- **Diseases**: Follow standard terminologies (as per MONDO)
- **GO Terms & Pathways**: Follow standard terminologies (e.g., as per Gene Ontology or Reactome)
- **Phenotypes**: First letter capitalized (e.g., `Calcified ovarian cyst`)
- **Proteins**: All uppercase (e.g., `TP53`)



