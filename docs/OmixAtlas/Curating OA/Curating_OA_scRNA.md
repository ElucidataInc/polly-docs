# Single Cell RNASeq Data

## Raw unfiltered counts


**I. Dataset Level Metadata**

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| abstract | abstract | Abstract | text | Abstract of the publication associated with the dataset. | Not available | Polly - Curated |
| author_cell_type | author_cell_type | Author Cell Type | text | Cell types as extracted from the publication/source | None | Polly - Curated |
| authors | authors | Authors | text | List of authors for the associated publication | None | Polly - Curated |
| curated_cell_capture_system | curated_cell_capture_system | Cell capture system | text | The cell capture system used in the protocol | None | Polly - Curated |
| curated_cell_line | kw_cell_line | Cell Line | text | Cell lines from which the samples were derived for this dataset | None | Polly - Curated |
| curated_chemistry_kit | curated_chemistry_kit | Chemistry kit | text | Chemistry/Reagent kit used in the protocol | None | Polly - Curated |
| curated_disease | disease | Disease | text | Disease associated with the dataset | Normal | Polly - Curated |
| curated_drug | kw_drug | Drug | text | Drugs administered to the samples belonging to the dataset | None | Polly - Curated |
| curated_library_protocol | curated_library_protocol | Library Protocol | text | Library preparation method used for this dataset | None | Polly - Curated |
| curated_modified_gene | curated_modified_gene | Gene | text | Gene(s) modified in the sample under study. | None | Polly - Curated |
| curated_organism | organism | Organism | text | The organism from which the samples were derived | Must be populated | Polly - Curated |
| curated_sampling_technique | curated_sampling_technique | Sampling Technique | text | The method/procedure used for collecting samples | None | Polly - Curated |
| curated_sequencer | curated_sequencer | Sequencer | text | The sequencing platform used in the protocol | None | Polly - Curated |
| curated_sequencing_technology | curated_sequencing_technology | Sequencing Technology | text | Sequencing method used in the protocol | None | Polly - Curated |
| curated_single_cell_chemistry | curated_single_cell_chemistry | Single cell Chemistry | text | The sequencing chemistry used for sequencing the single cell genome | None | Polly - Curated |
| curated_tissue | tissue | Tissue | text | Tissue from which the samples were derived | None | Polly - Curated |
| curated_treatment_name | curated_treatment_name | Treatment | text | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy | None | Polly - Curated |
| curated_treatment_response | curated_treatment_response | Treatment response | text | This field indicates the type/extent of response to the treatment | None | Polly - Curated |
| curated_treatment_type | curated_treatment_type | Treatment type | text | The type of treatment given to the sample | None | Polly - Curated |
| data_split_label | data_split_label | Data Split Label | text | Label on which data split was made. | None | Standard Identifier |
| data_type | kw_data_type | Data type | text | The type of biomolecular data captured (e.g. - Exome Sequencing, Transcriptomics) | Single cell | Standard Identifier |
| dataset_id | dataset_id | Dataset ID | text | Unique ID associated with every dataset | Must be populated | Standard Identifier |
| dataset_source | dataset_source | Source | text | Source from where the data was fetched | Must be populated | Standard Identifier |
| description | description | Description | text | Description of the dataset | Not available | Source Metadata |
| curated_donor_derived_dataset | curated_dataset_has_donor | Donor derived dataset | text | If the dataset is derived from a donor the value will be True else False | None | Polly - Curated |
| journal | journal | Journal | text | Journal in which the associated study was published | None | Polly - Curated |
| manually_curated_fields | manually_curated_fields | Manually curated fields | object | List of fields that have been manually curated for this dataset | text | Polly - Curated |
| overall_design | overall_design | Overall Design | text | Overall design of the experiment as given by the author | Not available | Source Metadata |
| publication_link | publication | Publication Link | text | Link to the publication associated with the dataset. If the associated publication information is not available, then this field provides the link to the data source providing more information on the dataset. | None | Polly - Curated |
| publication_title | publication_title | Publication Title | text | Title of the publication associated with this dataset | None | Polly - Curated |
| publication_year | year | Publication Year | integer | Year in which the dataset was published | Must be populated | Polly - Curated |
| pubmed_id | pubmed_ids | Pubmed ID | text | Pubmed unique identifier of the publication associated with the dataset | None | Polly - Curated |
| read_processing_pipeline | read_processing_pipeline | Processing pipeline | text | Software used in the original study to generate the raw counts matrix from sequencing data. | Not available | Standard Identifier |
| reference_annotation | reference_annotation | Reference annotation | text | Reference genome annotation used in the original study for counting/transcript quantification | Not available | Standard Identifier |
| reference_genome | reference_genome | Reference Genome | text | Reference genome build used for raw data processing at source to generate the counts matrix | Not available | Standard Identifier |
| source_link | source_link | Source Link | text | Link to data source from where the data was fetched | None | Polly - Curated |
| summary | summary | Summary | text | Summary of the experiment | Not available | Source Metadata |
| total_num_cells | total_num_cells | Cells | integer | Total number of cells in a dataset | integer | Polly - Curated |
| total_num_samples | total_num_samples | Samples | integer | Total number of samples in a dataset | integer | Polly - Curated |


**II Sample Level Metadata**

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| author_cell_type | kw_curated_raw_cell_type | Author cell type | text | Cell types as curated from the publication/source | none | Polly - Curated |
| cell_id | kw_column | Cell ID | text | Unique ID associated with every cell | Must be populated | Standard Identifier |
| curated_age_unit | curated_age_unit | Age unit | text | The age unit of the organism from which samples have been obtained. It is years for samples from humans and weeks for samples from mice. | none | Polly - Curated |
| curated_cell_line | kw_curated_cell_line | Cell line | text | Cell line from which the sample was derived | none | Polly - Curated |
| curated_developmental_stage_unit | curated_developmental_stage_unit | Developmental stage unit | text | Represents the unit used for representing the development stage | none | Polly - Curated |
| curated_developmental_stage_value | curated_developmental_stage_value | Developmental stage | text | Represents the stage of development/formation of the embryo of humans and mice | none | Polly - Curated |
| curated_disease | kw_curated_disease | Disease | text | Disease associated with the sample | normal | Polly - Curated |
| curated_disease_stage_system | curated_disease_stage_system | Disease stage type | text | This field represents the type of staging system indicating the disease stage | none | Polly - Curated |
| curated_disease_stage_value | curated_disease_stage_value | Disease stage | text | This field indicates the value of the disease stage | none | Polly - Curated |
| curated_donor_sample_type | curated_donor_sample_type | Donor sample type | text | The type of sample eg. Organoid, tissue, stem cell, etc. | none | Polly - Curated |
| curated_donor_type | curated_donor_type | Donor type | text | The type/clinical condition of the donor | none | Polly - Curated |
| curated_drug | kw_curated_drug | Drug | text | Drug administered in the sample | none | Polly - Curated |
| curated_gender | curated_gender | Gender | text | Gender of the organism from which the sample was derived | none | Polly - Curated |
| curated_genetic_mod_type | kw_curated_genetic_mod_type | Genetic Modification | text | Signifies the kind of genetic modification done on the sample. | none | Polly - Curated |
| curated_max_age | curated_max_age | Maximum age | text | The upper limit of the age range of the organism from which the samples have been obtained for the study. | none | Polly - Curated |
| curated_min_age | curated_min_age | Minimum age | text | The lower limit of the age range of the organism from which the samples have been obtained for the study. | none | Polly - Curated |
| curated_modified_gene | curated_gene_modified | Modified gene | text | Gene(s) modified in the sample under study. | none | Polly - Curated |
| curated_sampling_site | curated_sampling_site | Sampling site | text | The location/area from where the samples are collected. | none | Polly - Curated |
| curated_tissue | kw_curated_tissue | Tissue | text | Tissue from which the sample was taken | none | Polly - Curated |
| curated_treatment_name | curated_treatment_name | Treatment | text | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy | none | Polly - Curated |
| curated_treatment_response | curated_treatment_response | Treatment response | text | This field indicates the type/extent of response to the treatment | none | Polly - Curated |
| curated_treatment_type | curated_treatment_type | Treatment type | text | The type of treatment given to the sample | none | Polly - Curated |
| n_genes_by_counts | n_genes_by_counts | Genes per cell | float | Number of genes detected per cell | Must be populated | Standard Identifier |
| pct_counts_mt | pct_counts_mt | Mitochondrial proportion | float | Proportion of total counts for a cell which is mitochondrial | Must be populated | Standard Identifier |
| sample_id | sample | Sample ID | text | Unique ID associated with a sample | Must be populated | Standard Identifier |
| title | title | Title | text | Title of the sample | none | Source Metadata |
| total_counts | total_counts | Total counts | float | Absolute number of observed transcripts for a cell | Must be populated | Standard Identifier |


**III Feature Level Metadata**

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| total_counts | total_counts | Total counts | float | Sum of counts for feature over all cells | Must be populated | Polly - Curated |
| pct_dropout_by_counts | pct_dropout_by_counts | Percentage dropouts | float | Percentage of cells that gene does not appear in | Must be populated | Polly - Curated |
| n_cells | n_cells | Total cells | integer | Number of cells containing the gene | Must be populated | Polly - Curated |
| mean_counts | mean_counts | Mean counts | float | Average counts per cell for feature | Must be populated | Polly - Curated |
| gene_id | gene_id | Gene ID | text | NCBI/Entrez gene ID | none | Standard Identifier |
| feature_id | feature_id | Feature ID | text | ID of the feature (gene, protein, metabolite, etc.) being measured | Must be populated | Standard Identifier |
| ensembl_id | ensembl_id | Ensemble ID | text | Ensemble ID | none | Standard Identifier |


## Polly Processed Data
**I. Dataset Level Metadata**

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| abstract | abstract | Abstract | text | Abstract of the publication associated with the dataset. | Not available | Polly - Curated |
| author_cell_type | author_cell_type | Author Cell Type | text | Cell types as extracted from the publication/source | None | Polly - Curated |
| authors | authors | Authors | text | List of authors for the associated publication | None | Polly - Curated |
| curated_cell_capture_system | curated_cell_capture_system | Cell capture system | text | The cell capture system used in the protocol | None | Polly - Curated |
| curated_cell_line | kw_cell_line | Cell Line | text | Cell lines from which the samples were derived for this dataset | None | Polly - Curated |
| curated_chemistry_kit | curated_chemistry_kit | Chemistry kit | text | Chemistry/Reagent kit used in the protocol | None | Polly - Curated |
| curated_disease | disease | Disease | text | Disease associated with the dataset | Normal | Polly - Curated |
| curated_drug | kw_drug | Drug | text | Drugs administered to the samples belonging to the dataset | None | Polly - Curated |
| curated_library_protocol | curated_library_protocol | Library Protocol | text | Library preparation method used for this dataset | None | Polly - Curated |
| curated_modified_gene | curated_modified_gene | Gene | text | Gene(s) modified in the sample under study. | None | Polly - Curated |
| curated_organism | organism | Organism | text | The organism from which the samples were derived | Must be populated | Polly - Curated |
| curated_sampling_technique | curated_sampling_technique | Sampling Technique | text | The method/procedure used for collecting samples | None | Polly - Curated |
| curated_sequencer | curated_sequencer | Sequencer | text | The sequencing platform used in the protocol | None | Polly - Curated |
| curated_sequencing_technology | curated_sequencing_technology | Sequencing Technology | text | Sequencing method used in the protocol | None | Polly - Curated |
| curated_single_cell_chemistry | curated_single_cell_chemistry | Single cell Chemistry | text | The sequencing chemistry used for sequencing the single cell genome | None | Polly - Curated |
| curated_tissue | tissue | Tissue | text | Tissue from which the samples were derived | None | Polly - Curated |
| curated_treatment_name | curated_treatment_name | Treatment | text | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy | None | Polly - Curated |
| curated_treatment_response | curated_treatment_response | Treatment response | text | This field indicates the type/extent of response to the treatment | None | Polly - Curated |
| curated_treatment_type | curated_treatment_type | Treatment type | text | The type of treatment given to the sample | None | Polly - Curated |
| data_split_label | data_split_label | Data Split Label | text | Label on which data split was made. | None | Standard Identifier |
| data_type | kw_data_type | Data type | text | The type of biomolecular data captured (e.g. - Exome Sequencing, Transcriptomics) | Single cell | Standard Identifier |
| dataset_id | dataset_id | Dataset ID | text | Unique ID associated with every dataset | Must be populated | Standard Identifier |
| dataset_source | dataset_source | Source | text | Source from where the data was fetched | Must be populated | Standard Identifier |
| description | description | Description | text | Description of the dataset | Not available | Source Metadata |
| curated_donor_derived_dataset | curated_dataset_has_donor | Donor derived dataset | text | If the dataset is derived from a donor the value will be True else False | None | Polly - Curated |
| journal | journal | Journal | text | Journal in which the associated study was published | None | Polly - Curated |
| manually_curated_fields | manually_curated_fields | Manually curated fields | object | List of fields that have been manually curated for this dataset | text | Polly - Curated |
| overall_design | overall_design | Overall Design | text | Overall design of the experiment as given by the author | Not available | Source Metadata |
| publication_link | publication | Publication Link | text | Link to the publication associated with the dataset. If the associated publication information is not available, then this field provides the link to the data source providing more information on the dataset | None | Polly - Curated |
| publication_title | publication_title | Publication Title | text | Title of the publication associated with this dataset | None | Polly - Curated |
| publication_year | year | Publication Year | integer | Year in which the dataset was published | Must be populated | Polly - Curated |
| pubmed_id | pubmed_ids | Pubmed ID | text | Pubmed unique identifier of the publication associated with the dataset | None | Polly - Curated |
| read_processing_pipeline | read_processing_pipeline | Processing pipeline | text | Software used in the original study to generate the raw counts matrix from sequencing data. | Not available | Standard Identifier |
| reference_annotation | reference_annotation | Reference annotation | text | Reference genome annotation used in the original study for counting/transcript quantification | Not available | Standard Identifier |
| reference_genome | reference_genome | Reference Genome | text | Reference genome build used for raw data processing at source to generate the counts matrix | Not available | Standard Identifier |
| source_link | source_link | Source Link | text | Link to data source from where the data was fetched | None | Polly - Curated |
| summary | summary | Summary | text | Summary of the experiment | Not available | Source Metadata |
| total_num_cells | total_num_cells | Cells | integer | Total number of cells in a dataset | integer | Polly - Curated |
| total_num_samples | total_num_samples | Samples | integer | Total number of samples in a dataset | integer | Polly - Curated |
| batch_corrected | batch_corrected | Batch Corrected | boolean | Displays information whether the batch correction was performed ( False is batch_correction is not done ) | False | Polly - Curated |
| batch_correction_method | batch_correction_method | Batch Correction Method | text | Displays information on what batch correction was performed | none | Polly - Curated |
| curated_markers | curated_markers | Curated Markers | object | List of cell types: aggregated marker | Must be populated | Polly - Curated |
| embeddings_available | embeddings_available | Embeddings | text | List of embeddings available with the dataset | Must be populated | Polly - Curated |
| polly_curated_cell_type | polly_curated_cell_type | Polly Curated Cell Type | text | Types of cells present in the dataset, curated with standard ontology | Must be populated | Polly - Curated |
| polly_raw_cell_type | polly_raw_cell_type | Polly processed raw Cell type | text | Raw cell type associated with the cell | Must be populated | Polly - Curated |


**II Sample Level Metadata**

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| author_cell_type | kw_curated_raw_cell_type | Author cell type | text | Cell types as curated from the publication/source | none | Polly - Curated |
| cell_id | kw_column | Cell ID | text | Unique ID associated with every cell | Must be populated | Standard Identifier |
| curated_age_unit | curated_age_unit | Age unit | text | The age unit of the organism from which samples have been obtained. It is years for samples from humans and weeks for samples from mice. | none | Polly - Curated |
| curated_cell_line | kw_curated_cell_line | Cell line | text | Cell line from which the sample was derived | none | Polly - Curated |
| curated_developmental_stage_unit | curated_developmental_stage_unit | Developmental stage unit | text | Represents the unit used for representing the development stage | none | Polly - Curated |
| curated_developmental_stage_value | curated_developmental_stage_value | Developmental stage | text | Represents the stage of development/formation of the embryo of humans and mice | none | Polly - Curated |
| curated_disease | kw_curated_disease | Disease | text | Disease associated with the sample | normal | Polly - Curated |
| curated_disease_stage_system | curated_disease_stage_system | Disease stage type | text | This field represents the type of staging system indicating the disease stage | none | Polly - Curated |
| curated_disease_stage_value | curated_disease_stage_value | Disease stage | text | This field indicates the value of the disease stage | none | Polly - Curated |
| curated_donor_sample_type | curated_donor_sample_type | Donor sample type | text | The type of sample eg. Organoid, tissue, stem cell, etc. | none | Polly - Curated |
| curated_donor_type | curated_donor_type | Donor type | text | The type/clinical condition of the donor | none | Polly - Curated |
| curated_drug | kw_curated_drug | Drug | text | Drug administered in the sample | none | Polly - Curated |
| curated_gender | curated_gender | Gender | text | Gender of the organism from which the sample was derived | none | Polly - Curated |
| curated_genetic_mod_type | kw_curated_genetic_mod_type | Genetic Modification | text | Signifies the kind of genetic modification done on the sample. | none | Polly - Curated |
| curated_max_age | curated_max_age | Maximum age | text | The upper limit of the age range of the organism from which the samples have been obtained for the study. | none | Polly - Curated |
| curated_min_age | curated_min_age | Minimum age | text | The lower limit of the age range of the organism from which the samples have been obtained for the study. | none | Polly - Curated |
| curated_modified_gene | curated_gene_modified | Modified gene | text | Gene(s) modified in the sample under study. | none | Polly - Curated |
| curated_sampling_site | curated_sampling_site | Sampling site | text | The location/area from where the samples are collected. | none | Polly - Curated |
| curated_tissue | kw_curated_tissue | Tissue | text | Tissue from which the sample was taken | none | Polly - Curated |
| curated_treatment_name | curated_treatment_name | Treatment | text | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy | none | Polly - Curated |
| curated_treatment_response | curated_treatment_response | Treatment response | text | This field indicates the type/extent of response to the treatment | none | Polly - Curated |
| curated_treatment_type | curated_treatment_type | Treatment type | text | The type of treatment given to the sample | none | Polly - Curated |
| n_genes_by_counts | n_genes_by_counts | Genes per cell | float | Number of genes detected per cell | Must be populated | Standard Identifier |
| pct_counts_mt | pct_counts_mt | Mitochondrial proportion | float | Proportion of total counts for a cell which is mitochondrial | Must be populated | Standard Identifier |
| sample_id | sample | Sample ID | text | Unique ID associated with a sample | Must be populated | Standard Identifier |
| title | title | Title | text | Title of the sample | none | Source Metadata |
| total_counts | total_counts | Total counts | float | Absolute number of observed transcripts for a cell | Must be populated | Standard Identifier |
| pct_counts_ribo | pct_counts_ribo | Ribosomal percentage | float | Percentage counts per cell from ribosomal genes | Must be populated | Polly - Curated |
| pct_counts_in_top_20_genes | pct_counts_in_top_20_genes | Top 20 gene counts percentage | float | Percentage of total counts per cell contributed by top 20 most expressed genes | Must be populated | Polly - Curated |
| pct_counts_hb | pct_counts_hb | Hemo percentage | float | Percentage counts per cell from hemo genes | Must be populated | Polly - Curated |
| doublet_score | doublet_score | Doublet score | float | Predicted probability of a cell being doublet; output of a doublet prediction tool | Must be populated | Polly - Curated |
| clusters | clusters | Clusters | integer | Number associated with each cluster from clustering algorithm | Must be populated | Polly - Curated |
| polly_raw_cell_type | polly_raw_cell_type | Polly processed raw Cell type | text | Raw cell type associated with the cell | none | Polly - Curated |
| polly_curated_cell_type | polly_curated_cell_type | Polly processed standardized Cell type | text | Curated cell type associated with the cell | none | Polly - Curated |
| curated_cell_ontology_id | curated_cell_ontology_id | Cell ontology ID | text | Curated Cell Ontology ID | none | Polly - Curated |


**III Feature Level Metadata**

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| total_counts | total_counts | Total counts | float | Sum of counts for feature over all cells | Must be populated | Polly - Curated |
| pct_dropout_by_counts | pct_dropout_by_counts | Percentage dropouts | float | Percentage of cells that gene does not appear in | Must be populated | Polly - Curated |
| n_cells | n_cells | Total cells | integer | Number of cells containing the gene | Must be populated | Polly - Curated |
| mean_counts | mean_counts | Mean counts | float | Average counts per cell for feature | Must be populated | Polly - Curated |
| gene_id | gene_id | Gene ID | text | NCBI/Entrez gene ID | none | Standard Identifier |
| feature_id | feature_id | Feature ID | text | ID of the feature (gene, protein, metabolite etc.) being measured | Must be populated | Standard Identifier |
| ensembl_id | ensembl_id | Ensemble ID | text | Ensemble ID | none | Standard Identifier |
| means | means | Mean | float | mean value per gene | Must be populated | Polly - Curated |
| dispersions | dispersions | Dispersion | float | dispersions per gene | Must be populated | Polly - Curated |
| dispersions_norm | dispersions_norm | Normalized dispersion | float | normalized dispersions per gene | Must be populated | Polly - Curated |
| var | var | Variance | float | variance per gene | Must be populated | Polly - Curated |
| var_norm | var_norm | Variance dispersion | float | variance dispersions per gene | Must be populated | Polly - Curated |
| highly_variable | highly_variable | Highly variable gene | text | Whether a gene is highly variable | Must be populated | Polly - Curated |
| highly_variable_nbatches | highly_variable_nbatches | Highly variable in nbatches | integer | no. of batches/samples in which the gene is identified as highly variable | Must be populated | Polly - Curated |


## Author processed counts

 **I. Dataset Level Metadata**

| **Field** | **Description** | **Ontology** | **GUI Display Name** | **Polly-Python Display Name** |
| --- | --- | --- | --- | --- |
| Organism | This field represents the organism from which the samples originated. Organism labels already present in the source metadata are normalized using a normalization model. In case the organism labels are missing, related texts and abstracts are processed and further normalized to get the organism metadata label. | [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) | Organism | curated\_organism |
| Tissue | This field represents the tissue(s) from which the samples in the dataset are derived. Tissue labels already present in the source metadata are normalized using a normalization model. In cases where tissue labels are missing, related texts and abstracts are processed and further normalized to get the tissue metadata label. <br/> Tissue labels for datasets consists of all tissue names from which the samples are derived. Tissue labels are annotated for samples extracted from a healthy tissue or a diseased tissue. <br/> Key specifications for tissue metadata annotations are as follows: <ul><li> All labels are harmonized with Brenda Tissue Ontology </li><li> Dual Channel datasets where the author has studied two different cases in a single sample are not curated with tissue metadata </li><li> Datasets with numerical metadata are not annotated | [Brenda Tissue Ontology](https://www.ebi.ac.uk/ols/ontologies/bto) | Tissue | curated\_tissue |
| Drug | This field represents the drug(s) that have been used in the treatment of the samples or relate to the experiment in some other way. Drug labels already present in the source metadata are normalized using a normalization model. In cases where drug labels are missing, related texts and abstracts are processed and further normalized to get the drug metadata label. <br/> Drug labels are annotated for the following types of sample treatments: <ul><li> Chemical treatment - Treatment for activation/regulation/inhibition of a protein or nucleic acid. </li><li> Stimulation - Treatment to elicit an immune response. Eg. Cytokine stimulation (Interferon alpha, TNF-a), LPS </li><li> Drug treatment - Treatment with a substance used to treat an illness, relieve a symptom or modify a chemical process in the body for a specific purpose. Eg. 3,5-diethoxycarbonyl-1,4-dihydrocolidine (DDC) </li></ul> <br/> Drug labels are not annotated for the following types of sample treatments: <ul><li> Control/ Vehicle - The control group receives either no treatment, a standard treatment whose effect is already known, or a placebo. Mostly the control sample contains organic solvents </li><li> Treatment for genetic perturbation - Treatments used for clonal selection, inducible genetic perturbation. Eg. Puromycin, Tamoxifen, Doxycycline </li><li> Transfection - RNA molecules (miRNA, shRNA, siRNA etc) for genetic manipulation </li><li> Other treatments - Treatments such as chemotherapy, radiation therapy or with antibodies </li><li> Other - Culture media, supplements and detergents. Eg. DMEM, LB Media, Agar Media, Glucose, amino acids, fats, lipids, SDS </li></ul> <br/> **Note** : Any mention of the drug in the text is included as a drug label irrespective of whether it is being used in the experiment or not. | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Drug | curated\_drug |
| Disease | This field represents the disease(s) being studied in the experiment. Disease labels already present in the source metadata are normalized using a normalization model. In case the disease labels are missing, related texts and abstracts are processed and further normalized to get the disease metadata label. <br/> Disease labels are annotated when the samples have been collected from diseased tissue or organism. Examples of such cases are as follows: <ul><li> If the cell line has been extracted from a diseased organism/tissue and been immortalised for further study then it will be annotated for the disease. For example, if the samples are extracted from an osteosarcoma tumor then disease annotation for this dataset will be "osteosarcoma" </li><li> If a disease model has been created through genetic modification or has been artificially created in lab for further studies, then it will be annotated for that disease. For example, a genetically engineered mouse model where an immunodeficient or humanized mouse is implanted with a tissue from a patient's tumor. </li></ul> <br/> **Note** : In studies, where the cell lines are extracted from a healthy tissue/organism and then conditioned to induce disease, the disease label for such a dataset will be "normal" since the sample is not extracted from any diseased tissue or organism. <br/> Key specifications for disease metadata annotations are: <ul><li> All types of disease including viral, and bacterial infections, metabolic syndromes such as obesity, and diabetes, cancers such as lung neoplasm, breast neoplasm etc. are included as disease metadata labels. </li><li> Labels for disease mentioned for both full form and short form are annotated. For example, Acute Myeloid Leukaemia will be labelled for both AML, Myeloid Leukaemia, Leukaemia and Acute Myeloid Leukaemia </li><li> For a dataset, diseases for each sample in the dataset are annotated </li><li> Disease mentioned in the metadata of a dataset, irrespective of the study is labelled </li><li> Processes such as "carcinogenesis" or "tumorigenesis" are not annotated as diseases </li><li> Development abnormalities are not included </li><li> Dual Channel datasets where the author has studied two different cases in a single sample are not curated for disease metadata | [MeSH](https://www.ncbi.nlm.nih.gov/mesh/) | Disease | curated\_disease |
| Cell type | This field represents the cell type of the samples within the study. <br/> Cell type labels are annotated in cases where the authors have cultured a particular cell type either extracted from tissues or developmental organs or generated in the lab and then used it in further experiment. <br/> Key specifications for cell type metadata annotations are as follows: <ul><li> All labels are harmonized with Cell Ontology </li><li> Cell lines are not curated for cell types </li><li> For cell types with functional terms such as circulating, associated, derived, etc., only the cell type is annotated </li><li> For tissue-specific cell types, the tissue name along with the cell type is labelled. Eg - aortic endothelial cells, spinal motor neurons </li><li> Organism terms are not included in the cell type labels. Eg. mouse "HSPCs" </li><li> Abbreviated cell types with functional conditions are labelled as abbreviated terms. Eg. CTCs </li></ul> | [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl) | Cell Type | curated\_cell\_type |
| Cell line | This field represents the cell line from which the samples were extracted. List of the population of modified cells used for the study. Cell line labels already present in the source metadata are normalized using synonyms present in the cell line ontology we use. <br/> Cell line labels are annotated in cases where the authors have cultured a particular cell line or bought it from organizations such as ATCC and then used in further experiments. Eg. MDAMB-231, HEK-293 <br/> Key specifications for cell line metadata annotations are as follows: <ul><li> All labels are harmonized with The Cellosaurus Ontology </li><li> Dual Channel datasets where the author has studied two different cases in a single sample are not curated for cell line metadata </li><li> Datasets with numerical metadata are not annotated </li><ul> | [The Cellosaurus](https://web.expasy.org/cellosaurus/) | Cell Line | curated\_cell\_line |
| **Other metadata fields** | | |
| **Field** | **Description** | **GUI Display Name** | **Polly-Python Display Name** |
| Abstract | This field provides the abstract of the publication associated with the dataset. | NA | abstract |
| Year | This field provides the year in which the dataset or study is published. | Year | year |
| Gene | This field provides the gene(s) studied in the dataset. | Gene | curated_gene |
| Single cell chemistry | This field represents the sequencing method/platform used for sequencing the single cell genome. | Single cell chemsitry | curated_single_cell_chemistry |
| Sampling technique | This field represents the method/procedure used for collecting samples. | Sampling technique | curated_sampling_technique |
| Sample storage technique | This field represents the method/technique used for preservation or storage of sample for analysis. | Sampling storage technique | curated_storage_technique |
| Summary | This field provides a detailed summary of the publication (can be the abstract) or a summary of the experiment. | Summary (Available for datasets from GEO only) | summary |
| Overall design | This field provides information on the overall design of the experiment as given by the author. | Overall Design (Available for datasets from GEO only) | overall_design |
| Publication | This field provides the link to the publication associated with the dataset. If the associated publication information is not available, then this field provides the link to the data source providing more information on the dataset. | NA | publication |
| Source | This field provides the name of the source repository from where the dataset is fetched. | Source | dataset_source |
| Description | This field provides a brief description of the experiment or the study. | <ul><li> Title (view details page & card view) </li><li> Description (Table view) | description |
| Data Type | This field provides the type of biomolecular data represented/studied in the dataset. | NA | data_type |
| Dataset ID | This field provides the unique id for the dataset/study to represent a group of samples. | Dataset ID | dataset_id |
| Number of Cells | This field represents the number of cells/observations in the dataset. | Number of Cells | total_num_cells |
| Number of Samples | This field represents the total number of samples in a dataset. | Samples | total_num_samples |

#### II. Sample Level Metadata

| **Field** | **Description** | **Ontology** | **GUI Display Name** | **Polly-Python Display Name** |
| --- | --- | --- | --- | --- |
| **Ontology-driven Fields** |
| Tissue | This field represents the tissue(s) from which the samples originated. Tissue labels already present in the source metadata are normalized using a normalization model. In cases where tissue labels are missing, related texts and abstracts are processed and further normalized to get the tissue metadata label. <br/>Tissue labels are annotated for samples extracted from healthy or diseased tissue. All labels are harmonized with Brenda Tissue Ontology. | [Brenda Tissue Ontology](https://www.ebi.ac.uk/ols/ontologies/bto) | Tissue | curated\_tissue |
| Disease | At the sample level, this field represents the disease associated with a particular sample. Disease labels already present in the source metadata are normalized using a normalization model. In case the disease labels are missing, related texts and abstracts are processed and further normalized to get the disease metadata label. <br/> Disease labels are annotated for a sample when the samples have been collected from diseased tissue or organism. Examples of such cases are as follows: <ul><li> If the cell line has been extracted from a diseased organism/tissue and been immortalised for further study then it will be annotated for the disease. For example, if the samples under study are extracted from osteosarcoma tumors then disease annotation for this sample will be "osteosarcoma" </li><li>If a disease model has been created through genetic modification or has been artificially created in a lab for further studies, then it will be annotated for that disease. For example, a genetically engineered mouse model where an immunodeficient or humanized mouse is implanted with tissue from a patient's tumor. </li></ul> <br/> At the sample level, disease labels are annotated for the following sample type: <ul><li> Clinical- Samples extracted from diseased patients, tissue, cell lines etc. </li><li> GEM (Genetically Engineered Models) - Samples extracted from genetically engineered mouse models </li><li> Diet Induced - Samples extracted from diet induced mouse models </li><li> Xenograft - Samples extracted from patient-derived xenograft or cell line-derived xenograft mouse models </li><li> Infection - Samples extracted from infected organisms, tissues, cell lines or cultures. Example: Viral infection, Bacterial infection etc </li><li> Other - Any other type of sample in which disease has not been mentioned in metadata, but is not normal </li></ul> | [MeSH](https://www.ncbi.nlm.nih.gov/mesh/) | Disease | curated\_disease |
| Drug | This field represents the drugs that have been used in the treatment of a sample. Drug labels already present in the source metadata are normalized using a normalization model. In cases where drug labels are missing, related texts and abstracts are processed and and further normalized to get the drug metadata label. <br/> Drug labels are annotated for the following types of sample treatments: <ul><li> Chemical treatment - Treatment for activation/regulation/inhibition of a protein or nucleic acid. </li><li> Stimulation - Treatment to elicit an immune response. Eg. Cytokine stimulation (Interferon alpha, TNF-a), LPS </li><li> Drug treatment - Treatment with a substance used to treat an illness, relieve a symptom or modify a chemical process in the body for a specific purpose. Eg. 3,5-diethoxycarbonyl-1,4-dihydrocolidine (DDC) </li></ul> <br/> Drug labels are not annotated for the following types of sample treatments: <ul><li> Control/ Vehicle - The control group receives either no treatment, a standard treatment whose effect is already known, or a placebo. Mostly the control sample contains organic solvents </li><li> Treatment for genetic perturbation - Treatments used for clonal selection, inducible genetic perturbation. Eg. Puromycin, Tamoxifen, Doxycycline </li><li> Transfection - RNA molecules (miRNA, shRNA, siRNA etc) for genetic manipulation </li><li> Other treatments - Treatments such as chemotherapy, radiation therapy or with antibodies </li><li> Other - Culture media, supplements and detergents. Eg. DMEM, LB Media, Agar Media, Glucose, amino acids, fats, lipids, SDS </li></ul> <br/> Note: Any mention of the drug in the text is included as a drug label irrespective of whether it is being used in the experiment or not. | [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Drug | curated\_drug |
| Cell line | This field represents the cell line from which the sample was derived. Cell line labels already present in the source metadata are normalized using synonyms present in the cell line ontology we use. The cell line field is curated for a sample if the authors have cultured a particular cell line or bought it from organisations such as ATCC and then used in the further experiment. The names of the cell lines are harmonized by the cellosaurus ontology. | [The Cellosaurus](https://web.expasy.org/cellosaurus/) | Cell Line | curated\_cell\_line |
| Cell Type | This field represents the cell type of the sample. Cell type labels are annotated where the authors have cultured a particular cell type either extracted from tissues or developmental organs or generated in the lab and then used it in the further experiment.The cell type field provides the closest cell type name as per the Cell Ontology. This cell type label can be either source derived or by manual cell type annotation. | [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl) | Cell Type | curated\_cell\_type |
| **Other Metadata Fields** |
| **Field** | **Description** | **GUI Display Name** | **Polly-Python Display Name** |
| Gene | Gene of interest in the sample | Gene | curated\_gene |
| Genetic Modification | This fields represents the kind of genetic modification done on the sample. | Genetic Modification | curated\_genetic\_modification\_type |
| Modified Genes | This fields represents the gene(s) modified in the sample under study. | NA | curated\_gene\_modified |
| Donor Type | This field represents the type/clinical condition of the donor | Donor Type | curated\_donor\_type |
| Donor Sample Type | This field represents the location/area from where the tumor samples are collected | Donor Sample Type | curated\_donor\_sample\_type |
| Gender | This fields represents the gender of the organism from which the sample was derived | Gender | curated\_gender |
| Minimum age | This fields provides the lower limit of the age range of the organism from which the samples have been obtained for the study. | Minimum age | curated\_min\_age |
| Maximum age | This fields provides the upper limit of the age range of the organism from which the samples have been obtained for the study. | Maximum age | curated\_max\_age |
| Age unit | This fields provides the age unit of the organism from which samples have been obtained. It is years for samples from humans and weeks for samples from mice. | Age unit | curated\_age\_unit |
| Sampling site | This fields provides the location/area from where the tumor samples are collected. | Sampling site | curated\_sampling\_site |
| Treatment | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy | Treatment Name | curated\_treatment\_name |
| Treatment type | The type of treatment given to the sample | Treatment Type | curated\_treatment\_type |
| Response to Treatment | This field indicates the type/extent of response on the treatment | Treatment Response | curated\_treatment\_response |
| Author Cell Type | This field represents the author cell type as mentioned in the publication associated with the dataset. | Author cell type | curated\_raw\_cell\_type |
| Marker present | This field represents the gene name/names that are differentially expressed in a cluster based on which the cell type of the cluster is annotated. | Marker Present | curated\_marker\_present |
| Marker absent | This field represents the gene name/names that are absent in a cluster based on which the cell type of the cluster is annotated. | Marker Absent | curated\_marker\_absent |
| Cell ontology ID | This field represents the unique ID for the cell type according to the Cell Ontology. | NA | curated\_cell\_ontology\_id |
| Clusters | This field provides the cluster number to which each cell belongs after subjecting the dataset to the clustering process (Using Leiden or other algorithms). | Cell Type Cluster | clusters |
| UMI Counts | This field represents the Unique Molecular Identifier count per cell. It represents an absolute number of observed transcripts. The number should be higher than 500 in a cell. | UMI Counts | umi\_counts |
| Sample ID | This field represents the unique ID of the sample. | Sample ID | sample\_id |
| Cell ID | This field represents the unique ID associated with every cell | NA | cell\_id |
| Gene Counts | This field represents the number of genes detected per cell (defined as genes with cpm \> 1); cpm=counts per million | Gene Counts | gene\_counts |
| Mitochondrial count | This field represents the percentage of mitochondrial counts in total counts for a cell. | NA | percent\_mito |
| Title | This field represents the title of the sample, representing the type/genotype/origin/experimental condition of the sample | NA | title |

#### III. Feature Level Metadata

| **Field** | **Description** | **Polly-Python Display Name** |
| --- | --- | --- |
| Feature ID | This field represents the ID of the feature (gene, metabolite, protein etc) being measured. | feature_id |
| Highly variable gene | This field indicates whether the gene is highly variable. <br>For highly variable genes - True; otherwise -False. | highly_variable |
| Number of cells | This field represents the number of cells which are containing the gene. | n_cells |

### Cell Type Annotation for single cell datasets: Manual Curation Process

Cell-type labels are assigned at the cell cluster level based on expression signatures using ontology or controlled vocabularies. For datasets, where cell-type annotations are not available from the source (mainly GEO datasets), we manually curate the cell-type information based on the differential marker expression for clusters. In cases where cell type annotations are already available in datasets at source, datasets are not manually re-curated.

![Process Flow](../img/OmixAtlas-Images/1_C.png) <center>**Figure 1.** Process Flow</center>

#### **I) Manual identification of the cell types and markers from Publications** -

Internal curators determine if a particular dataset can be curated for cell type by going through the publications associated with the dataset. In publications, the information on cell type and the corresponding marker is present either in the figures (UMAP, T-SNE plots), text or supplementary files. If this information is not present, then such datasets are marked as 'Not Curatable'. The following types of studies fall under the category of 'Not curatable' datasets.

- Single Cell Type Study - The whole study was done on only one cell type.
- Lineages - The study included the lineage of one cell type. Example - T helper cells, T memory cells etc.
- Cell Cycle Studies - In this study the differential markers were studied for the G1, S, G2 and M phases of the cell cycle for a particular cell type.
- Methods - Different analysis methods were studied
- Publications having \< 2000 cells - These publications used methods which could not be reproduced. Therefore, these were not curated.
- Publication not available
- Marker information not available in the publication
  - Marker Cell Type Info Absent
  - Marker Info Absent
- Time Point
- Cell lines
- Transitional Cells
- Cell Type Info Absent
- Embryonic Development
- Organoids
- Others

#### II) Metadata addition- Annotation of clusters

The process of cell-type cluster annotation for curatable datasets is based on the general scRNASeq general[workflow](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) using the Scanpy library with steps as shown below in the figure:
                             
![Homepage](../img/OmixAtlas-Images/2_C.png) <center>**Figure 2.** Workflow</center>

UMAP/tSNE plots are generated as a result of single-cell raw count processing. By visualization of clusters with UMAP/t-SNE plots, cell type cluster annotation is done.

##### 1. Cluster annotation with raw cell type (cell type terminology used in publication):

Based on the marker expression value for each cluster, the cell type is annotated to the cluster. This annotation is added as a field named curated\_raw\_cell\_type. The raw cell type cluster annotation is compared with cell type annotation from the publication such as:

- All or most cell types are annotated
- UMAP is structurally similar to the one in the publication
- Relative proportions of cell types are matching
- Relative positions of cell types are matching
- Ontological terms and marker information are added for the cell type

##### 2. Ontological terms and marker information:

1. Cell type ontology + ontology ID: Cell-type annotation corresponding to Cell Ontology. This is given as the field named: curated\_cell\_type
2. Marker information: Gene name/names that are differentially expressed in the cluster
  - curated\_marker\_present: Gene name/names that are differentially expressed in the cluster
  - curated\_marker\_absent: Gene name/names that are absent in the cluster
  - 

### Cell Type Visualization

Users can see the relative distribution of cell numbers across different cell types when viewing the details of a single-cell RNA-seq dataset before selecting a dataset for analysis.

![Homepage](../img/OmixAtlas-Images/OA_view_5.png) 

**NOTE** : For manual cell type curation, datasets should be available on Polly.
                             

