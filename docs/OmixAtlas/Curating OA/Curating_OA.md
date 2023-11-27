
## What is Curation?

Curation is a process of making biological data organized, annotated, and standardized to make it a FAIR (Findable, Accessible, Interoperable, and Reusable) resource.

- **Findable** - to assign a persistent identifier that would make data easy to find by both humans and machines
- **Accessible** - to be able to retrieve the data by their identifier using a standard protocol
- **Interoperable** - to use standardized terms and have references to other data and be machine actionable
- **Reusable** - to sufficiently describe data for both computers and humans to be able to understand the data

Using a structured FAIR data/database can fast-track the drug development process as it make it easier for researchers to test their first-hand hypotheses.


## How is Data Curated on OA?

Curated datasets stored on Polly are a result of a multi-step process that ensures the overall quality and reusability of the data. In summary, all data sets on Polly go through the following two-step process:-

- Data Engineering - This includes transforming data to fit a proprietary data schema that is uniform across several datatypes. The transformation streamlines data in one consistent schema and allows users to query multiple data types on a single data infrastructure.
- Harmonizing Metadata - This involves tagging each sample and dataset with a metadata harmonized with uniform ontology.


## Why is Data Curation Important?

Manually curating vast volumes of unstructured or semi-structured biomedical data for drug development can be expensive, cumbersome, time-consuming, and resource-intensive. Automating the curation process would ensure the following -

- Ensuring high data quality.
- Enabling data reusability for various types of analysis.
- Facilitating the rapid curation of recently added datasets.

## Curated Fields on Polly

On Polly, all datasets and the corresponding samples are available with six standard metadata fields which are curated using Polly's proprietary NLP-based Polly-BERT models and harmonized with specific biomedical ontologies. Such fields are Organism, Disease, Tissue, Cell Line, Cell Type and Drug. In addition to these fields, other metadata fields are also available which are extracted from the source metadata. Curated metadata on Polly is available at 3 levels:

- Dataset-level metadata - General information about the experiment, subject and transformations, for eg, organism, experiment type, and disease under study.
- Sample-level metadata - Captures information for each sample. Eg. Drug, tissue
- Feature-level metadata: Provides molecular/gene information that is consistent across samples.
