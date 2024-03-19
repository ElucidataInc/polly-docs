# Single Cell RNASeq Data

## Raw unfiltered counts


### I. Dataset Level Metadata

| Field name (Polly-python)   | Original Name (Data)   | Display Name (GUI)  | Type   | Description | Default values | Metadata Category |
|---------|-------|-----------|--------|-------|----------------|-----------|
| abstract                 | abstract                 | Abstract                  | text   | Abstract of the publication associated with the dataset.                                                        | Not available  | Source Metadata   |
| author_cell_type         | author_cell_type         | Author Cell Type          | text   | Cell types as extracted from the publication/source                                                             | None           | Polly - Curated   |
| curated_author_cell_type| curated_author_cell_type| Curated Author Cell Type  | text   | Cell types as extracted from the publication/source and standardized using cell ontology                        | None           | Polly - Curated   |
| authors                  | authors                  | Authors                   | text   | List of authors for the associated publication                                                                    | None           | Polly - Curated   |
| curated_cell_capture_system| curated_cell_capture_system| Cell capture system     | text   | The cell capture system used in the protocol                                                                    | None           | Polly - Curated   |
| curated_cell_line        | kw_cell_line             | Cell Line                 | text   | Cell lines from which the samples were derived for this dataset'                                                 | None           | Polly - Curated   |
| curated_chemistry_kit    | curated_chemistry_kit    | Chemistry kit             | text   | Chemistry/Reagent kit used in the protocol                                                                     | None           | Polly - Curated   |
| curated_disease          | disease                  | Disease                   | text   | Disease associated with the dataset                                                                             | Normal         | Polly - Curated   |
| curated_drug             | kw_drug                  | Drug                      | text   | Drugs administered to the samples belonging to the dataset                                                      | None           | Polly - Curated   |
| curated_library_protocol | curated_library_protocol | Library Protocol          | text   | Library preparation method used for this dataset                                                               | None           | Polly - Curated   |
| curated_modified_gene    | curated_modified_gene    | Gene                      | text   | Gene(s) modified in the sample under study.                                                                     | None           | Polly - Curated   |
| curated_organism         | organism                 | Organism                  | text   | The organism from which the samples were derived                                                               | Must be populated | Polly - Curated   |
| curated_sampling_technique| curated_sampling_technique| Sampling Technique       | text   | The method/procedure used for collecting samples                                                               | None           | Polly - Curated   |
| curated_sequencer        | curated_sequencer        | Sequencer                 | text   | The sequencing platform used in the protocol                                                                   | None           | Polly - Curated   |
| curated_sequencing_technology| curated_sequencing_technology| Sequencing Technology| text   | Sequencing method used in the protocol'                                                                        | None           | Polly - Curated   |
| curated_single_cell_chemistry| curated_single_cell_chemistry| Single cell Chemistry| text   | The sequencing chemistry used for sequencing the single cell genome                                            | None           | Polly - Curated   |
| curated_tissue           | tissue                   | Tissue                    | text   | Tissue from which the samples were derived'                                                                    | None           | Polly - Curated   |
| curated_treatment_name   | curated_treatment_name   | Treatment                 | text   | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy                              | None           | Polly - Curated   |
| curated_treatment_response| curated_treatment_response| Treatment response        | text   | This field indicates the type/extent of response to the treatment                                              | None           | Polly - Curated   |
| curated_treatment_type   | curated_treatment_type   | Treatment type            | text   | The type of treatment given to the sample'                                                                     | None           | Polly - Curated   |
| data_split_label         | data_split_label         | Data Split Label          | Array  | Label on which data split was made.                                                                             | None           | Standard Identifier |
| data_type                | kw_curated_data_type     | Data type                 | text   | The type of biomolecular data captured (e.g. - Exome Sequencing, Transcriptomics)                                | Single cell    | Standard Identifier |
| dataset_id               | dataset_id               | Dataset ID                | text   | Unique ID associated with every dataset'                                                                       | Must be populated | Standard Identifier |
| dataset_source           | dataset_source           | Source                    | text   | Source from where the data was fetched'                                                                        | Must be populated | Standard Identifier |
| description              | description              | Description               | text   | Description of the dataset'                                                                                    | Not available  | Source Metadata   |
| curated_donor_derived_dataset| curated_dataset_has_donor| Donor derived dataset  | text   | If the dataset is derived from a donor the value will be True else False                                        | None           | Polly - Curated   |
| journal                  | journal                  | Journal                   | text   | Journal in which the associated study was published                                                            | None           | Polly - Curated   |
| manually_curated_fields  | manually_curated_fields  | Manually curated fields   | object | List of fields that have been manually curated for this dataset                                                 | text           | Polly - Curated   |
| overall_design           | overall_design           | Overall Design            | text   | Overall design of the experiment as given by the author                                                         | Not available  | Source Metadata   |
| publication_link         | publication              | Publication Link          | text   | Link to the publication associated with the dataset. If the associated publication information is not available, then this field provides the link to the data source providing more information on the dataset.| None           | Polly - Curated   |
| publication_title        | publication_title        | Publication Title         | text   | Title of the publication associated with this dataset                                                           | None           | Polly - Curated   |
| publication_year         | year                     | Publication Year          | integer| Year in which the dataset was published                                                                         | Must be populated | Polly - Curated   |
| pubmed_id                | pubmed_ids               | Pubmed ID                 | text   | Pubmed unique identifier of the publication associated with the dataset                                         | None           | Polly - Curated   |
| read_processing_pipeline | read_processing_pipeline | Processing pipeline       | text   | Software used in the original study to generate the raw counts matrix from sequencing data.                     | Not available  | Standard Identifier |
| reference_annotation     | reference_annotation     | Reference annotation      | text   | Reference genome annotation used in the original study for counting/transcript quantification                   | Not available  | Standard Identifier |
| reference_genome         | reference_genome         | Reference Genome          | text   | Reference genome build used for raw data processing at source to generate the counts matrix                     | Not available  | Standard Identifier |
| source_link              | source_link              | Source Link               | text   | Link to data source from where the data was fetched                                                            | None           | Polly - Curated   |
| summary                  | summary  | Summary | text   | Summary of the experiment | Not available  | Source Metadata   |
| total_num_cells          | total_num_cells          | Cells | integer | Total number of cells in a dataset | integer | Polly - Curated |
| total_num_samples | total_num_samples | Samples | integer | Total number of samples in a dataset | integer|   |



### II Sample Level Metadata

| Field name (Polly-python)   | Original Name (Data)   | Display Name (GUI)  | Type   | Description | Default values | Metadata Category |
|-------------------------|-------|------------------------|--------|-----------|----------------|-------------------|
| author_cell_type        | kw_curated_raw_cell_type | Author cell type       | Array  | Cell types as curated from the publication/source                                                      | none           | Polly - Curated   |
| curated_author_cell_type | curated_author_cell_type | Curated Author Cell Type | Array  | Cell types as extracted from the publication/source and standardized using cell ontology             | None           | Polly - Curated   |
| cell_id                 | kw_column                | Cell ID                | text   | Unique ID associated with every cell'                                                                  | Must be populated | Standard Identifier |
| curated_age_unit        | curated_age_unit         | Age unit               | text   | The age unit of the organism from which samples have been obtained. It is years for samples from humans and weeks for samples from mice. | none           | Polly - Curated   |
| curated_cell_line       | kw_curated_cell_line     | Cell line              | text   | Cell line from which the sample was derived                                                           | none           | Polly - Curated   |
| curated_developmental_stage_unit | curated_developmental_stage_unit | Developmental stage unit | text   | Represents the unit used for representing the development stage                                        | none           | Polly - Curated   |
| curated_developmental_stage_value | curated_developmental_stage_value | Developmental stage | text   | Represents the stage of development/formation of the embryo of humans and mice                         | none           | Polly - Curated   |
| curated_disease         | kw_curated_disease       | Disease                | text   | Disease associated with the sample'                                                                   | normal         | Polly - Curated   |
| curated_disease_stage_system | curated_disease_stage_system | Disease stage type | text   | This field represents the type of staging system indicating the disease stage'                         | none           | Polly - Curated   |
| curated_disease_stage_value | curated_disease_stage_value | Disease stage         | text   | This field indicates the value of the disease stage                                                    | none           | Polly - Curated   |
| curated_donor_sample_type | curated_donor_sample_type | Donor sample type     | text   | The type of sample eg. Organoid, tissue, stem cell, etc.                                               | none           | Polly - Curated   |
| curated_donor_type      | curated_donor_type       | Donor type             | text   | The type/clinical condition of the donor                                                               | none           | Polly - Curated   |
| curated_drug            | kw_curated_drug          | Drug                   | text   | Drug administered in the sample                                                                        | none           | Polly - Curated   |
| curated_gender          | curated_gender           | Gender                 | text   | Gender of the organism from which the sample was derived                                               | none           | Polly - Curated   |
| curated_genetic_mod_type | kw_curated_genetic_mod_type | Genetic Modification | text   | Signifies the kind of genetic modification done on the sample.                                          | none           | Polly - Curated   |
| curated_max_age         | curated_max_age          | Maximum age            | text   | The upper limit of the age range of the organism from which the samples have been obtained for the study. | none           | Polly - Curated   |
| curated_min_age         | curated_min_age          | Minimum age            | text   | The lower limit of the age range of the organism from which the samples have been obtained for the study. | none           | Polly - Curated   |
| curated_modified_gene   | curated_gene_modified    | Modified gene          | text   | Gene(s) modified in the sample under study.                                                             | none           | Polly - Curated   |
| curated_sampling_site   | curated_sampling_site    | Sampling site          | text   | The location/area from where the samples are collected.                                                 | none           | Polly - Curated   |
| curated_tissue          | kw_curated_tissue        | Tissue                 | text   | Tissue from which the sample was taken                                                                  | none           | Polly - Curated   |
| curated_treatment_name  | curated_treatment_name   | Treatment              | text   | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy                      | none           | Polly - Curated   |
| curated_treatment_response | curated_treatment_response | Treatment response  | text   | This field indicates the type/extent of response to the treatment                                        | none           | Polly - Curated   |
| curated_treatment_type  | curated_treatment_type   | Treatment type         | text   | The type of treatment given to the sample'                                                              | none           | Polly - Curated   |
| n_genes_by_counts       | n_genes_by_counts        | Genes per cell         | float  | Number of genes detected per cell                                                                      | Must be populated | Standard Identifier |
| pct_counts_mt           | pct_counts_mt            | Mitochondrial proportion | float | Proportion of total counts for a cell which is mitochondrial                                             | Must be populated | Standard Identifier |
| sample_id               | sample                   | Sample ID              | text   | Unique ID associated with a sample                                                                     | Must be populated | Standard Identifier |
| title                   | title                    | Title                  | text   | Title of the sample | none  | Source Metadata   |
| total_counts | total_counts | Total counts | float  | Absolute number of observed transcripts for a cell | Must be populated | Standard Identifier |


### III Feature Level Metadata

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

### I. Dataset Level Metadata

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| abstract                | abstract                 | Abstract                | text    | Abstract of the publication associated with the dataset.                                                   | Not available     | Polly - Curated   |
| author_cell_type        | author_cell_type        | Author Cell Type        | text    | Cell types as extracted from the publication/source                                                        | None              | Polly - Curated   |
| curated_author_cell_type | curated_author_cell_type | Curated Author Cell Type | text    | Cell types as extracted from the publication/source and standardized using cell ontology                  | None              | Polly - Curated   |
| authors                 | authors                  | Authors                 | text    | List of authors for the associated publication                                                              | None              | Polly - Curated   |
| curated_cell_capture_system | curated_cell_capture_system | Cell capture system  | text    | The cell capture system used in the protocol                                                                | None              | Polly - Curated   |
| curated_cell_line       | kw_curated_cell_line     | Cell Line               | text    | Cell lines from which the samples were derived for this dataset                                             | None              | Polly - Curated   |
| curated_chemistry_kit   | curated_chemistry_kit    | Chemistry kit           | text    | Chemistry/Reagent kit used in the protocol                                                                 | None              | Polly - Curated   |
| curated_disease         | disease                  | Disease                 | text    | Disease associated with the dataset                                                                         | Normal            | Polly - Curated   |
| curated_drug            | kw_curated_drug          | Drug                    | text    | Drugs administered to the samples belonging to the dataset                                                  | None              | Polly - Curated   |
| curated_library_protocol | curated_library_protocol | Library Protocol        | text    | Library preparation method used for this dataset                                                            | None              | Polly - Curated   |
| curated_modified_gene   | curated_modified_gene    | Gene                    | text    | Gene(s) modified in the sample under study.                                                                 | None              | Polly - Curated   |
| curated_organism        | organism                 | Organism                | text    | The organism from which the samples were derived                                                            | Must be populated | Polly - Curated   |
| curated_sampling_technique | curated_sampling_technique | Sampling Technique   | text    | The method/procedure used for collecting samples                                                            | None              | Polly - Curated   |
| curated_sequencer       | curated_sequencer        | Sequencer               | text    | The sequencing platform used in the protocol                                                                | None              | Polly - Curated   |
| curated_sequencing_technology | curated_sequencing_technology | Sequencing Technology | text    | Sequencing method used in the protocol'                                                                     | None              | Polly - Curated   |
| curated_single_cell_chemistry | curated_single_cell_chemistry | Single cell Chemistry | text    | The sequencing chemistry used for sequencing the single cell genome                                         | None              | Polly - Curated   |
| curated_tissue          | tissue                   | Tissue                  | text    | Tissue from which the samples were derived'                                                                 | None              | Polly - Curated   |
| curated_treatment_name  | curated_treatment_name   | Treatment               | text    | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy                           | None              | Polly - Curated   |
| curated_treatment_response | curated_treatment_response | Treatment response    | text    | This field indicates the type/extent of response to the treatment                                            | None              | Polly - Curated   |
| curated_treatment_type  | curated_treatment_type   | Treatment type          | text    | The type of treatment given to the sample'                                                                  | None              | Polly - Curated   |
| data_split_label        | data_split_label         | Data Split Label        | text    | Label on which data split was made.                                                                         | None              | Standard Identifier |
| data_type               | kw_curated_data_type     | Data type               | text    | The type of biomolecular data captured (e.g. - Exome Sequencing, Transcriptomics)                           | Single cell       | Standard Identifier |
| dataset_id              | dataset_id               | Dataset ID              | text    | Unique ID associated with every dataset'                                                                    | Must be populated | Standard Identifier |
| dataset_source          | dataset_source           | Source                  | text    | Source from where the data was fetched'                                                                     | Must be populated | Standard Identifier |
| description             | description              | Description             | text    | Description of the dataset'                                                                                 | Not available     | Source Metadata   |
| curated_donor_derived_dataset | curated_dataset_has_donor | Donor derived dataset | text    | If the dataset is derived from a donor the value will be True else False                                     | None              | Polly - Curated   |
| journal                 | journal                  | Journal                 | text    | Journal in which the associated study was published                                                         | None              | Polly - Curated   |
| manually_curated_fields | manually_curated_fields | Manually curated fields | object  | List of fields that have been manually curated for this dataset                                             | text              | Polly - Curated   |
| overall_design          | overall_design           | Overall Design          | text    | Overall design of the experiment as given by the author                                                     | Not available     | Source Metadata   |
| publication_link        | publication              | Publication Link        | text    | Link to the publication associated with the dataset. If the associated publication information is not available, then this field provides the link to the data source providing more information on the dataset. | None | Polly - Curated |
| publication_title       | publication_title        | Publication Title       | text    | Title of the publication associated with this dataset                                                        | None              | Polly - Curated   |
| publication_year        | year                     | Publication Year        | integer | Year in which the dataset was published                                                                     | Must be populated | Polly - Curated   |
| pubmed_id               | pubmed_ids               | Pubmed ID               | text    | Pubmed unique identifier of the publication associated with the dataset                                      | None              | Polly - Curated   |
| read_processing_pipeline | read_processing_pipeline | Processing pipeline     | text    | Software used in the original study to generate the raw counts matrix from sequencing data.                 | Not available     | Standard Identifier |
| reference_annotation    | reference_annotation     | Reference annotation    | text    | Reference genome annotation used in the original study for counting/transcript quantification               | Not available     | Standard Identifier |
| reference_genome        | reference_genome         | Reference Genome        | text    | Reference genome build used for raw data processing at source to generate the counts matrix                  | Not available     | Standard Identifier |
| source_link             | source_link              | Source Link             | text    | Link to data source from where the data was fetched                                                         | None              | Polly - Curated   |
| summary                 | summary                  | Summary                 | text    | Summary of the experiment                                                                                   | Not available     | Source Metadata   |
| total_num_cells         | total_num_cells          | Cells                   | integer | Total number of cells in a dataset'                                                                         | integer           | Polly - Curated   |
| total_num_samples       | total_num_samples        | Samples                 | integer | Total number of samples in a dataset'                                                                       | integer           | Polly - Curated   |
| batch_corrected         | batch_corrected          | Batch Corrected         | boolean | Displays information whether the batch correction was performed ( False is batch_correction is not done )    | False             | Polly - Curated   |
| batch_correction_method | batch_correction_method  | Batch Correction Method| text    | Displays information on what batch correction was performed                                                 | none              | Polly - Curated   |
| curated_markers         | curated_markers          | Curated Markers         | object  | List of cell types: aggregated marker                                                                       | Must be populated | Polly - Curated   |
| embeddings_available    | embeddings_available     | Embeddings              | text    | List of embeddings available with the dataset                                                               | Must be populated | Polly - Curated   |
| polly_curated_cell_type | polly_curated_cell_type  | Polly Curated Cell Type | text    | Types of cells present in the dataset, curated with standard ontology                                        | None              | Polly - Curated   |
| polly_raw_cell_type     | polly_raw_cell_type      | Polly processed raw Cell type | text | Raw cell type associated with the cell                                                                      | None              | Polly - Curated   |


### II Sample Level Metadata

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| author_cell_type       | kw_curated_raw_cell_type | Author cell type          | text   | Cell types as curated from the publication/source                                                | none              | Polly - Curated   |
| curated_author_cell_type | curated_author_cell_type | Curated Author Cell Type | text   | Cell types as extracted from the publication/source and standardized using cell ontology       | None              | Polly - Curated   |
| cell_id                | kw_curated_column      | Cell ID                     | text   | Unique ID associated with every cell                                                           | Must be populated | Standard Identifier |
| curated_age_unit       | curated_age_unit       | Age unit                    | text   | The age unit of the organism from which samples have been obtained.                             | none              | Polly - Curated   |
| curated_cell_line      | kw_curated_cell_line   | Cell line                   | text   | Cell line from which the sample was derived                                                    | none              | Polly - Curated   |
| curated_developmental_stage_unit | curated_developmental_stage_unit | Developmental stage unit | text   | Represents the unit used for representing the development stage                                | none              | Polly - Curated   |
| curated_developmental_stage_value | curated_developmental_stage_value | Developmental stage    | text   | Represents the stage of development/formation of the embryo of humans and mice                  | none              | Polly - Curated   |
| curated_disease        | kw_curated_disease     | Disease                     | text   | Disease associated with the sample                                                             | normal            | Polly - Curated   |
| curated_disease_stage_system | curated_disease_stage_system | Disease stage type    | text   | This field represents the type of staging system indicating the disease stage                   | none              | Polly - Curated   |
| curated_disease_stage_value | curated_disease_stage_value | Disease stage            | text   | This field indicates the value of the disease stage                                             | none              | Polly - Curated   |
| curated_donor_sample_type | curated_donor_sample_type | Donor sample type       | text   | The type of sample eg. Organoid, tissue, stem cell, etc.                                        | none              | Polly - Curated   |
| curated_donor_type     | curated_donor_type     | Donor type                  | text   | The type/clinical condition of the donor                                                        | none              | Polly - Curated   |
| curated_drug           | kw_curated_drug        | Drug                        | text   | Drug administered in the sample                                                                  | none              | Polly - Curated   |
| curated_gender         | curated_gender         | Gender                      | text   | Gender of the organism from which the sample was derived                                          | none              | Polly - Curated   |
| curated_genetic_mod_type | kw_curated_genetic_mod_type | Genetic Modification  | text   | Signifies the kind of genetic modification done on the sample.                                    | none              | Polly - Curated   |
| curated_max_age        | curated_max_age        | Maximum age                 | text   | The upper limit of the age range of the organism from which the samples have been obtained for the study. | none         | Polly - Curated   |
| curated_min_age        | curated_min_age        | Minimum age                 | text   | The lower limit of the age range of the organism from which the samples have been obtained for the study. | none         | Polly - Curated   |
| curated_modified_gene  | curated_gene_modified  | Modified gene               | text   | Gene(s) modified in the sample under study.                                                      | none              | Polly - Curated   |
| curated_sampling_site  | curated_sampling_site  | Sampling site               | text   | The location/area from where the samples are collected.                                           | none              | Polly - Curated   |
| curated_tissue         | kw_curated_tissue      | Tissue                      | text   | Tissue from which the sample was taken                                                           | none              | Polly - Curated   |
| curated_treatment_name | curated_treatment_name | Treatment                   | text   | Name of the treatment given to the samples i.e. name of the chemical/drug/therapy               | none              | Polly - Curated   |
| curated_treatment_response | curated_treatment_response | Treatment response    | text   | This field indicates the type/extent of response to the treatment                                  | none              | Polly - Curated   |
| curated_treatment_type | curated_treatment_type | Treatment type               | text   | The type of treatment given to the sample                                                         | none              | Polly - Curated   |
| n_genes_by_counts      | n_genes_by_counts      | Genes per cell              | float  | Number of genes detected per cell                                                                | Must be populated | Standard Identifier |
| pct_counts_mt          | pct_counts_mt          | Mitochondrial proportion    | float  | Proportion of total counts for a cell which is mitochondrial                                      | Must be populated | Standard Identifier |
| sample_id              | sample                  | Sample ID                   | text   | Unique ID associated with a sample                                                                | Must be populated | Standard Identifier |
| title                  | title                   | Title                       | text   | Title of the sample                                                                             | none              | Source Metadata   |
| total_counts           | total_counts           | Total counts                | float  | Absolute number of observed transcripts for a cell                                                 | Must be populated | Standard Identifier |
| pct_counts_ribo        | pct_counts_ribo        | Ribosomal percentage        | float  | Percentage counts per cell from ribosomal genes                                                   | Must be populated | Polly - Curated   |
| pct_counts_in_top_20_genes | pct_counts_in_top_20_genes | Top 20 gene counts percentage | float | Percentage of total counts per cell contributed by top 20 most expressed genes                     | Must be populated | Polly - Curated   |
| pct_counts_hb          | pct_counts_hb          | Hemo percentage             | float  | Percentage counts per cell from hemo genes                                                        | Must be populated | Polly - Curated   |
| doublet_score          | doublet_score          | Doublet score               | float  | Predicted probability of a cell being doublet; output of a doublet prediction tool                  | Must be populated | Polly - Curated   |
| clusters               | clusters               | Clusters                    | integer| Number associated with each cluster from clustering algorithm                                      | Must be populated | Polly - Curated   |
| polly_raw_cell_type    | polly_raw_cell_type    | Polly processed raw Cell type | text  | Raw cell type associated with the cell                                                            | none              | Polly - Curated   |
| polly_curated_cell_type | polly_curated_cell_type | Polly processed standardized Cell type | text | Curated cell type associated with the cell                                                        | none              | Polly - Curated   |
| polly_curated_cell_ontology_id | curated_cell_ontology_id | Cell ontology ID       | text   | Curated Cell Ontology ID                                                                         | none              | Polly - Curated   |

### III Feature Level Metadata

| Field name(Polly-Python) | Original Name(Data) | Display Name (GUI) | Type | Description | Default values | Metadata Category |
|----|------|----|------|---------|-------|-------|
| total_counts             | total_counts          | Total counts              | float   | Sum of counts for feature over all cells              | Must be populated | Polly - Curated   |
| pct_dropout_by_counts    | pct_dropout_by_counts | Percentage dropouts       | float   | Percentage of cells that gene does not appear in      | Must be populated | Polly - Curated   |
| n_cells                  | n_cells               | Total cells               | integer | Number of cells containing the gene                   | Must be populated | Polly - Curated   |
| mean_counts              | mean_counts           | Mean counts               | float   | Average counts per cell for feature                   | Must be populated | Polly - Curated   |
| gene_id                  | gene_id               | Gene ID                   | text    | NCBI/Entrez gene ID                                   | none              | Standard Identifier |
| feature_id               | feature_id            | Feature ID                | text    | ID of the feature (gene, protein, metabolite etc.) being measured | Must be populated | Standard Identifier |
| ensembl_id               | ensembl_id            | Ensemble ID               | text    | Ensemble ID                                           | none              | Standard Identifier |
| means                    | means                 | Mean                      | float   | mean value per gene                                   | Must be populated | Polly - Curated   |
| highly_variable         | highly_variable      | Highly variable gene      | text    | Whether a gene is highly variable                     | Must be populated | Polly - Curated   |
| highly_variable_nbatches | highly_variable_nbatches | Highly variable in nbatches | integer | no. of batches/samples in which the gene is identified as highly variable | Must be populated | Polly - Curated   |



### Cell Type Annotation for single cell datasets: Manual Curation Process

Cell-type labels are assigned at the cell cluster level based on expression signatures using ontology or controlled vocabularies. For datasets, where cell-type annotations are not available from the source (mainly GEO datasets), we manually curate the cell-type information based on the differential marker expression for clusters. In cases where cell type annotations are already available in datasets at source, datasets are not manually re-curated.


**I) Manual identification of the cell types and markers from Publications** 

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

**II) Metadata addition- Annotation of clusters**

The process of cell-type cluster annotation for curatable datasets is based on the general scRNASeq general[workflow](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) using the Scanpy library with steps as shown below in the figure:
                            

UMAP/tSNE plots are generated as a result of single-cell raw count processing. By visualization of clusters with UMAP/t-SNE plots, cell type cluster annotation is done.

1. Cluster annotation with raw cell type (cell type terminology used in publication):

Based on the marker expression value for each cluster, the cell type is annotated to the cluster. This annotation is added as a field named curated\_raw\_cell\_type. The raw cell type cluster annotation is compared with cell type annotation from the publication such as:

- All or most cell types are annotated
- UMAP is structurally similar to the one in the publication
- Relative proportions of cell types are matching
- Relative positions of cell types are matching
- Ontological terms and marker information are added for the cell type

2. Ontological terms and marker information:
   - Cell type ontology + ontology ID: Cell-type annotation corresponding to Cell Ontology. This is given as the field named: curated\_cell\_type
   - Marker information: Gene name/names that are differentially expressed in the cluster
      - curated\_marker\_present: Gene name/names that are differentially expressed in the cluster
      - curated\_marker\_absent: Gene name/names that are absent in the cluster
    

### Cell Type Visualization

Users can see the relative distribution of cell numbers across different cell types when viewing the details of a single-cell RNA-seq dataset before selecting a dataset for analysis.

**NOTE** : For manual cell type curation, datasets should be available on Polly.
                             

