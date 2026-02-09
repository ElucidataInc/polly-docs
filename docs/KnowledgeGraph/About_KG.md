# Introduction

## Overview

**Polly Knowledge Graph** **(Polly KG)** is a cutting-edge knowledge graph designed to help researchers efficiently access and analyze internal and published data. It enables users to generate hypotheses and gain a comprehensive understanding of data related to specific targets, indications, and drug classes. By structuring relationships between biological entities, Polly KG eliminates siloed knowledge bases and streamlines data retrieval, empowering researchers to explore multi-modal evidence with  ease.

### Key Challenges Addressed

- **Manual Data Sifting**: Traditionally, researchers manually search through internal data and external literature to validate hypotheses, which is time-consuming and inefficient.
- **Siloed Knowledge Bases**: The lack of structured relationships between biological entities (e.g., genes, diseases, drugs) leads to fragmented information.
- **Limited Holistic Insights**: Researchers struggle to gather and integrate multi-modal evidence to support their biological questions.

## Key Components

**Polly Co-Scientist**

Polly Co-Scientist is an AI-powered research assistant that converts natural language questions into Cypher queries to interact with the Knowledge Graph. It allows users to query, explore relationships, and extract insights without requiring knowledge of Cypher or programming.

**Polly KG Explorer**

Polly KG Explorer is an interactive application within the Polly ecosystem designed to explore and analyze the Knowledge Graph. It enables users to visualize relationships between biological entities, apply filters, extract subgraphs, and generate data-driven hypotheses from harmonized biomedical data.

**Polly Python**

Polly Python is the programmatic interface for interacting with Polly, enabling users to access data, query Knowledge Graph resources, manage datasets, and run analyses using Python. It provides a flexible environment for automation, advanced querying, and integration into bioinformatics workflows.

### Terminology in the Polly KG Application

This section provides definitions of key terms used in the Polly KG Explorer application. Understanding these terms will help users navigate and utilize the tool more effectively.

### General Terms

- **Knowledge Graph (KG):** A structured representation of relationships between biological entities such as genes, proteins, diseases, and drugs. It enables data integration and the discovery of new insights.
- **Node:** A fundamental unit in the knowledge graph representing a biological entity (e.g., a gene, protein, disease, or drug).
- **Edge:** A connection between two nodes that represents a relationship between them (e.g., a gene encodes a protein).
- **Hops:** The number of steps (edges) between two nodes in the graph. A 1-hop relationship indicates that nodes are directly connected, while a 2-hop relationship indicates a connection through an intermediate node.
- **Searched Node:** The primary node that the user searches for in the application. This node serves as the central point for exploring relationships within the knowledge graph.
- **Immediate Neighbors:** Nodes that share a direct edge with the searched node. These nodes represent entities that have a direct relationship with the searched entity (e.g., a drug directly linked to a target protein).

### Application-Specific Terms

- **Node Selection Panel:** The section of the interface where users can select node categories and specific values for visualization.
- **Knowledge Graph Visualization Area:** The main workspace where selected nodes and their relationships are displayed as an interactive graph.
- **Polly Co-Scientist:** An AI-powered research assistant that transforms natural language queries into Cypher commands to interact with Polly’s Knowledge Graph, without requiring users to write Cypher code.

---

## Support

For additional assistance or questions related to the Knowledge Graph Query API, please contact the Elucidata support team.
