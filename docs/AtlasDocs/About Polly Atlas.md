# What is an Polly Atlas?

Polly Atlas is a user-defined collection of tables that integrates spreadsheet-style data with relational database structures. Each table is organized around a predefined set of metadata fields or clinical annotations, and provides unified access to linked datasets, QC reports, and other associated files or folders in a unified table.

Polly Atlas is designed to simplify how biomedical data is organized, accessed, and analyzed. It provides a structured, AI-ready environment where multi-modal datasets are harmonized, linked, and made easily queryable. Within this ecosystem, Polly Python enables programmatic access to Atlases, allowing users to search, retrieve, and work with curated datasets seamlessly in their analytical workflows.

Polly Atlas provides two complementary interfaces: **GUI** and **Code** designed to serve different user needs while enabling seamless interaction with the platform.

**GUI Interface**

This Polly Atlas provides an interface for basic exploration of the data, schema management and access management. Designed for intuitive, no-code interaction:

- **Data Exploration & Querying:** Filter, search, and explore datasets with ease  
- **Data Management:** Perform basic data operations through UI workflows  
- **Access Management:** Configure roles and permissions for secure data access  

**Code Interface**

Built for advanced users requiring flexibility and control via [Polly Python](https://docs.polly.elucidata.io/polly-python/Atlas/Atlas.html):

- **Data Management:** Programmatic handling of datasets at scale  
- **SQL Querying:** Advanced querying for complex analysis and transformations
- **Schema Management:** Define and manage table structures and metadata  


## Key Advantages

- **Eliminates Manual File Handling:** No need to browse S3/FTP folders or download files to explore data; everything is accessible within structured tables.  
- **Instant Data Discovery:** Search and filter datasets using column-level filters instead of scrolling through long file lists.  
- **Built-in Querying:** Perform complex filtering and data retrieval using intuitive UI filters or SQL queries.  
- **Relational Data Linking:** Seamlessly connect datasets, samples, and omics data through structured relationships; no more siloed files.  
- **Unified Data Access:** Access all file types (GCT, h5ad, STAR, Salmon outputs, reports) from a single interface. This platform provides unified data access to both primary and accessory files.  
- **Instant Data Snapshot:** Table overview, column descriptions, and ERD diagrams provide instant context without opening files externally.  
- **One-Click Data Export:** Download full/filtered tables instantly as CSV for downstream analysis.  
- **Full Audit Trail & Traceability:** Track all changes with user, action, and timestamp; ensuring reproducibility and compliance.  
- **Role-Based Access Control:** Granular permissions (Manager → Consumer) ensure secure and controlled data access.  
- **Improved User Experience:** Replaces fragmented workflows (download → script → filter) with a streamlined process: **Explore → Filter → Query → Download**







