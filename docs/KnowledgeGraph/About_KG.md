# Introduction

## Overview

**Polly Knowledge Graph** **(Polly KG)** is a cutting-edge knowledge graph designed to help researchers efficiently access and analyze internal and published data. It enables users to generate hypotheses and gain a comprehensive understanding of data related to specific targets, indications, and drug classes. By structuring relationships between biological entities, Polly KG eliminates siloed knowledge bases and streamlines data retrieval, empowering researchers to explore multi-modal evidence with  ease.

### Key Challenges Addressed

- **Manual Data Sifting**: Traditionally, researchers manually search through internal data and external literature to validate hypotheses, which is time-consuming and inefficient.
- **Siloed Knowledge Bases**: The lack of structured relationships between biological entities (e.g., genes, diseases, drugs) leads to fragmented information.
- **Limited Holistic Insights**: Researchers struggle to gather and integrate multi-modal evidence to support their biological questions.


## The Polly KG Explorer Application

The **Polly KG Explorer Application** is an interactive tool and a part of the **Polly KG ecosystem** on Polly. It is designed to help researchers explore and analyze the Polly Knowledge Base, which consists of harmonized biomedical data. It enables users to uncover relationships between biological entities and generate data-driven hypotheses.

### What Can Users Do with the Polly KG Explorer Application?

- **Explore Relationships**: Select specific genes, proteins, diseases, or drugs and visualize their connections within a structured knowledge graph.
- **Filter Views**: Use interactive filtering to refine queries and focus on relevant data.
- **Querying using Natural Language**: Query the Knowledge Graph using simple English, no need for Cypher or Python expertise. (Polly Co-Scientist)

### Additional Features of the Polly KG Explorer Application

- **Generate and Download knowledge graphs**: Extract specific portions of the KG, visualize them dynamically.
- **Navigate the KG with Ease**: Zoom, reset views, and switch between 2D and 3D visualizations for better insights.
- **Save & Share Insights**: Export graph images and download structured relationship data for further analysis.

### Who Should Use This Application?

- **Pharma R&D Teams**: Investigating potential drug targets and indications.
- **Multi-Omics Researchers**: Integrating and exploring diverse biological datasets to gain insights.
- **Bioinformatics Teams**: Conducting literature mining and entity linking to support discovery.

By offering a seamless, interactive experience, this application enables users to navigate complex biological relationships, extract meaningful insights, and accelerate hypothesis generation using the Polly KG Explorer framework.

---

### Terminology in the Polly KG Explorer Application

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
