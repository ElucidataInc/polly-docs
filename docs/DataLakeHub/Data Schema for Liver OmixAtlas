# Data Schema for Liver OmixAtlas

All types of data on Polly goes through Polly Curation process and are stored as tabular files or H5AD files (for single cell data) as described here. All data on Liver Omix Atlas have the metadata curated and harmonized using controlled vocabularies. There is data from various sources in the Liver OmixAtlas as follows :

1. GEO
2. TCGA
3. GTEx
4. MetaboLights
5. Metabolomics Workbench
6. Human Protein Atlas
7. CCLE
8. Depmap
9. LINCS
10. CPTAC

Here is the list of different data types available in the Liver OmixAtlas : 

1. Transcriptomics
2. Single Cell
3. Mutation
4. Metabolomics
5. Proteomics
6. Drug Screens
7. Gene Dependency
8. Gene Effect
9. Methylation
10. miRNA

Some metadata fields are common across all data types and sources whereas others are data type or source specific. The structure of data and metadata for each source and data type is described below:

## Common Metadata Fields Across Data Types

All the common metadata fields are curated in a consistent manner regardless of the source of the data. All the data in the Omix Atlas can be queried using these curated fields.

### 1.1 Dataset Level Fields 

All the datasets in the liver omix atlas has been annotated for unique id, description, organism, disease, tissue, experimental conditions, source, publication and curated for datatype, drug, cell line, cell type. These annotations helps in defining and mapping the characteristics of each data. These values has been standardized across all the datasets and are aligned with the FAIR guidelines. The dataset level mapping helps in narrowing down a dataset of interest.  

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| dataset_id  | String  | A unique id for dataset/study/project to represent a group of samples | ST000915_AN001489, GSE100155_GPL6884, LIHC_Proteomics_TCGA |
| description | String  | Brief text description providing details of the samples/experiment    | MiRNA profile of liver cell lines, Liver hepatocellular carcinoma methylation data |
| kw_data_type| String  | Identifies the type of measurement present in the dataset.  For example Transcriptomics data_type contains gene expression values, metabolomics contains intensity measurement from mass spec instruments etc. | Transcriptomics,Drug screens,Mutation |
| organism    | String/List | The organism over which the experiment for the dataset was conducted | Homo sapiens,Mus musculus |
| disease     | List | The name of the studied disease and associated effects | Fatty Liver, Inflammation, Obesity, Diabetes Mellitus |
| tissue      | List | The tissue on which the experiment was conducted for the said dataset | Liver, Adipose Tissue, Kidney, Brain | 
| kw_drug     | List | A database for drugs, chemical entity, natural and synthetic compounds which has been used as a treatment during the experiment in different samples | Troglitazone, Berberine, Clofibrate, Rosiglitazone |Studies with no drug perturbations will have this field as "None" |
| kw_cell_line | List | List of population of modified cells used for the study in quest | Hep-G2, MCF-7, HeLa, A549, TFK-1 | All TCGA datasets have this field as "None" |
| kw_cell_type | List | A list of differentitated cells which can be identified at morphological, structural and physiological level | Liver cells, Kidney Cells, Brain Cells |
| dataset_source | String | The name of the repository from where data has been originally deposited | GEO, TCGA, LINCS |
| publication | String | If the dataset has an associated publication, this field contains a link to the publication; in other cases, it contains a link to the data source providing more information regarding the dataset | 
| total_num_samples | Int | The total number of samples present in the dataset | 22, 84, 267 | 
| total_num_cells | Int | Number of cells in the experiment | 100 | This curated field is available only for Single cell studies | 


### 1.2 Sample Level Fields 

Sample level annotations directly defines the biological characteristics of each sample. All the samples in the liver omix atlas has been curated for disease, cell line, drug, cell type, genetic modifications, modified gene and tissue. The datasets can be queried using these fields as well. 

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| kw_column | String | A unique id for a sample | GSM173533, GSM4454526, HEP3B217_LIVER |
| kw_doc_id | String | Path to the dataset file | discover-prod-datalake-v1@@@liver_atlas@@data@@Transcriptomics@@GSE10409_GPL1355.gct |
| kw_curated_disease | String | Name of the disease condition for that particular sample | Normal, Liver Cirrhosis |
| kw_curated_cell_line | String | Population of modified cells used for the sample | SMMC-7721, Hep-G2, none | All TCGA samples have this field as "none" |
| kw_curated_drug | String | A database for drugs, chemical entity, natural and synthetic compounds which has been used as a treatment in the sample | Troglitazone, Berberine, Clofibrate, Rosiglitazone, | All the samples without drug perturbation will have this field as "none" |
| kw_curated_cell_type | String | Differentiated cells which can be identified at morphological, structural and physiological level | hepatocyte, endothelial cell of lymphatic vessel |
| kw_curated_genetic_mod_type | String | The type of genetic modification done on the sample  | wildtype, knockout, knockin, knockdown |
| kw_curated_modified_gene | String | A gene or list of genes modified in the sample | TP53, COL18A1, PRKAA2 | All the wildtype samples will have this field as "none" |
| kw_curated_tissue | String | A tissue or a list of tissues from which the sample has been obtained | Liver, Kidney, Brain, Colon |

### 1.3 Feature Level Fields 

The datasets has also been annotated for features of each sample which provides description about the gene being studied, unique identifier, path to the doc.

### 1.3.1 All Datatypes Except Single Cell  

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| kw_index | List | Symbol for the molecule (gene, metabolite, protein etc) studied in the experiment | Epas1, Tpgs1,Cxcl12 |
| kw_column.kw_column | List | Unique ID for each sample | GSM3034529, GSM2928029 |
| kw_column.kw_expression | List | Feature intensity of a sample. For genomics this will be expression value and for proteomics, lipidomics this will be metabolite intensity. | 0.7943000197410583, 0.3294999897480011, 0.6687999963760376 |
| kw_doc_id | List | Path to the dataset file. | discover-prod-datalake-v1@@@liver_atlas@@data@@Transcriptomics@@GSE10409_GPL1355.gct |

### 1.3.2 Single Cell

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| feature_name | List | Symbol for the gene studied in the experiment | Ntm |
| celltype | List | Unique ID for each sample | geneexp_cluster_7, geneexp_cluster_16, geneexp_cluster_4 |
| Value | List | Feature intensity of a sample | 0.000243683270913317, 0.000549465175768581 |
| Dataset | List | Path to the dataset | discover-prod-datalake-v1@@@liver_atlas@@data@@SingleCell@@GSE124395_GPL16791.h5ad |

### 1.4 Specific metadata fields from various sources

All the source specific fields which are mentioned in the following section are not curated by Polly and are present as they are in the source. Data in the Atlas can be queried using these fields as well. These fields may not be present in all the data on the source and hence as a result may not be present for all the data on Liver Omix Atlas as well.

### 1.4.1 CCLE (Sample level fields) 

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| ccle_id | List | This ID helps in mapping the cell line in CCLE database | HEP3B217_LIVER, SNU886_LIVER |
| depmapid | List | This is a unique ID which helps in finding the data in DepMap(Dependency Map) portal | ACH-000625, ACH-000739 |
| histology | List | This column provides the relevant information about the type of disease being studied | carcinoma, adenocarcinoma |
| name | List | Unique name for cell line | HEL9217, MCF7, 253JBV, SNU-878 |
| gender | List | Gender identity of the patient from whom the cell line has been obtained | Male, Female |
| age | List | Age of the patient from whom the cell line has been obtained | 8, 52, 56,43,28 |
| year | List | Year in which the sample was deposited | 2015, 2018 |

### 1.4.2 DepMap (Sample level fields) 

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| ccle_name | List | This list provides the detail about cell line and associated disease | HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE, LS513_LARGE_INTESTINE, MCF7_BREAST |
| cas9_activity | List | Percentage value for the efficiency of crispr cas 9 nuclease activity in the mammalian cell system | 65.2, 76.9, 92, 86.6, 52.4, 65.2 |
| sex | List | Gender identity of the patient from whom the cell line has been obtained | Male, Female |
| primary_disease | List | Original disease for which the cell line is model system | Kidney Cancer, Liver Cancer, Colon/Colorectal Cancer, Skin Cancer |
| subtype | List | A dictionary for smaller classes of cancer that a cancer can be grouped on the basis of physiological characteristics of cancer cell line | Melanoma, Adenocarcinoma, Rhabdomyosarcoma, Glioblastoma |
| age | List | Age of the patient from whom the cell line has been obtained | 8, 52, 56,43,28 |
| primary_tissue | List | Source tissue of cell line | Liver, Colon, Breast, Kidney |

### 1.4.3 LINCS (Sample level fields) 

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| pert_time | List | Depicts the numerical value of time of treatment of drugs during the experiment | 6, 12, 24, 48, 1 |
| pert_time_unit | List | Depicts the unit of time of treatment of drugs during the experiment | h, min, sec, week |
| pert_type | List | Column describing the experimental conditions eg: control/treated/untreated | ctl_vehicle, trt_cp, ctl_untrt |
| cell_id | List | Unique LINCS cell id for each cell line which acts as internal identifier in the LINCS database | A375.311, CL34 |
| pert_iname | List | Describes the condition of samples treated/untreated | UnTrt, Trt, DMSO |
| pert_dose | List | This column denotes the numerical value for the dose of chemical used during the experiment | 2, 25, 6 |
| pert_dose_unit | List | Provides information about the unit of the dose of the drug | uL, mL, mg |
| curated_is_control | List | Depicts information the sample is a control or a perturbation for a particular experiment | 1 |
| curated_cohort_id | List | The cohort of the sample | 2 |
| curated_cohort_name | List | Name of the group of samples | trt_cp - Gsk-429286a, ctl_vehicle - none |

### 1.4.4 MetaboLights (Sample level fields) 

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- |
| sample_source_name | List | Provides information about original condition/place/organization/laboratory from where the sample has been acquired. | BIIE Treatment Control, NPY Treatment, Moin Saleem, Bristol UK, IMTEK , S0_1_4_LAq, T5_2_9_MAq, plant |
| assay_ms_assay_name | List | Unique id which helps in mapping the samples origin/date of indexing/type of ms assay performed. | 0018_LC_20180917_sample_87, 0018_LC_20180917_sample_103, 0018_LC_20180917_sample_107 |
| subject id | List | Unique ID for sample factor value | DDO142, DDO93, DDO233 |
| factors_gender | List | Gender identity of the patient | Male, Female |
| factors_genotype | List | Contains information  | wild type genotype, Hi-MYC |
| factors_strain | List | Name of the strain of mice used during the experiment. | C57BL/6, B6 |
| factors_smoking status | List | Describes the smoking condition of the donor. | Never Smoker, Smoker |
| factors_cohort | List | Group of samples at a particular stage/in a particular experimental condition | validation, 1, 2 |
| factors_cell_type | List | Population of cells in which has been studied during the experiment. | control lung fibroblast, IPF lung fibroblast, Human peripheral blood mononuclear cells, non small cell lung cancer cell (NSCLC) |
| factors_plasma | List | Depicts the experimental condition of the samples which has been treated with plasma in a specific disease condition. |Non-thermal plasma sham, biological rep3, Time point 4 h after Non-thermal plasma treated 1 min, biological rep1, Cells Incubated with ACS-pre Patient Plasma, Cells Incubated with Normal Patient Plasma |
| factors_drug | List | Drugs used as a treatment. | WCB001_Mock_1d_RNA_4, DER, Tamoxifen, Deoxynivalenol, Solvent, Solvent+Spike |
| factors_injury | List | Provides information about type of injury of the donor. | alveolar cell injury, respiratory failure, AEC sham control, AEC mechanical injury 0 hr SW_N14, AEC cyclic stretch 8 hr SW_N25 |
| additional sample data_height cm | List | Height of the donor | 180.34, 164.084, 177.8 |
| additional sample data_weight kg | List | Weight of the donor | 102.965468, 80.73944186, 90.31024087 |
| additional sample data_tumor_type | List | Type of cancer donor is suffering from. | adenocarcinoma, lung adenocarcinoma, non-small cell lung cancer |
| factors_obesity | List | Type of obesity related disease donor is suffering from.  | insulin-sensitive (HOMA-IR<3) obese individual, liver dysfunction in obesity |

### 1.4.5 Metabolomics Workbench (Sample level fields)

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- | 
| sample_id | List | Unique Identity for each sample obtained from Metabolomics Workbench | NASH001, NASH013, NASH029, T26-6, T30-2 |
| Subject ID | List | Unique ID for the subject | SU0004318, SU0004310 |
| Factors.Diagnosis | List | List of disease which was diagnosed in the donor | Normal, Cirrhosis, Steatosis |
| Additional sample data.BMI | List | Body Mass Index of the patient | 43.3, 27, 34 |
| Additional sample data.AGE | List | Age of the patient | 45,72,83,34 |

### 1.4.6 GEO - Single Cell (Sample level fields)

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- | 
| umi_counts | Int | Unique Molecular Identifier count per cell. It represents absolute number of observed transcripts. The number should be higher than 500 in a cell. | 500-1000 |
| gene_counts | Int | Number of reads that align to each gene using different programs.  | 74.5446 fpkm, 30.6890 fpkm |
| sample | String | Unique Identifier for the sample.  | GSM2787809, GSM3396732 |
| characteristics_ch1 | List | Depicts the information about tissue/developmental stage/sample type/strain. | tissue: Pancreatic islets, developmental stage: adult, sample type: Single Cell |
| cell_type | List | This column lists the population of cells used during the experiment. | Dendritic cells, Basophils, CD4+ Cells, Progenitors, mesenchymal cells
| age | List | Age of the patient from whom the cell line has been obtained | 8, 52, 56,43,28 |
| donor | List | The organism from which the cells has been extracted  | Homo sapiens, Mus musculus |
| location | List | Origin place of the donor  | New York, USA, Bristol, UK |
| genes_detected | List | number of genes detected per sample (defined as genes with cpm > 1); cpm=counts per million | 
| donor_organism_sex | List | Gender identity of the donor organism | Male, Female |
| cell_subtype | List | Population of cells which has been studied during the experiment. | CD45+, GR1-, SSClow, CD11c+, MHCII+, CD11b+, CD24+, differentiated non lung cancer stem like cells (OSK-A549-SN) |
| cell type clust | List | Cluster of cell type which has been studied during the experiment. | C01_CD8-LEF1, C05_CD8-GZMK |
| gender | List | Gender identity of the donor organism | Male, Female |

### 1.4.7 GEO - Transcriptomics (Sample level fields)

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- | 
| geo_accession | List | Unique ID for each sample | GSM2699628, GSM2699632 |
| source_name_ch1 | List | Name of the source from which the sample has been obtained | liver cancer cells |
| curated_is_control | List | Depicts information the sample is a control or a perturbation for a particular experiment. | 1,2 |
| curated_cohort_id | List | The cohort of the sample. For one cohort the values will remain same. | 0,1 |
| curated_cohort_name | List | Name of the group of samples providing information about experimental condition/tissue/cell line/treatment. | Noodle diet; WAT_N-group; WAT_Noodle diet, hepatocellular carcinoma; liver cancer cells; hepg2_0Ã‚Âµg/ml_H |
| title | List | Description about the type/genotype/origin/experimental condition of the sample | hepg2_0Âµg/ml_H1 |

### 1.4.8 GTEX (Sample level fields)

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- | 
| SAMPID | List | Unique ID for each sample | GTEX-ZTPG-1426-SM-51MT3, GTEX-ZPU1-0826-SM-57WG2 |
| SMPTHNTS | List | Description about the treatment condition of the tissue | 2 pieces, diffuse macro and microvesucular steatosis, nodular regenerative |
| SMTS | List | Source Tissue | Liver |
| SMNABTCH | List | Batch ID | BP-43375, BP-43529 |
| DTHHRDY | List | List the cause of death of the patient | Fast death of natural causes, Ventilator Case, Intermediate death |

### 1.4.9 TCGA (Sample level fields)

| Field    | Type   |Description                     |Example values              | Exception            |
| :---------- | :---------- |:----------------------------------- | :----------------------------------- | :----------------------------------- | 
| sample_id | List | Unique ID for each sample | TCGA-2Y-A9GU-01A, TCGA-2V-A95S-01A |
| sample_type | List | The type of disease from which the sample has been extracted | Primary Tumor |
| primary_diagnosis | List | Type and cause of disease the patient has suffered from | Hepatocellular carcinoma-- NOS |
| primary_site | List | Primary site of tumor | Liver and intrahepatic bile ducts |
| disease_type | List | Type of disease the donor patient suffered from. | Adenomas and Adenocarcinomas |
| vital_status | List | Status of the donor patient. | Alive, Dead |
| bcr_patient_barcode | List | First four fields of the barcode. | TCGA-2Y-A9GU-01A, TCGA-2V-A95S-01A |
| barcode | List | Unique indexed identifier for each donor patient.  | TCGA-2Y-A9GU-01A-11R-A38B-07,  TCGA-2V-A95S-01A-11R-A37K-07 |
| patient_id | List | Unique ID of the donor. Contains first three fields of barcode.  | TCGA-2Y-A9GU, TCGA-2V-A95S |
| gender | List | Gender of the donor patient.  | Male, Female |
| age_at_diagnosis | List | Age of the patient during diagnosis. | 20187, 21318, 28387 |
| age_at_index | List | Age of the patient while indexing the data. | 55 |
| race | List | Racial classification of the donor patient. | white, caucasian, hispanic, latino | 
| Ethnicity | List | Cultural origin of the donor patient. | not hispanic or latino, African American, Caucasian American |
| primary site | List | Primary site of the tumor. | Liver and intrahepatic bile ducts |
| Ajcc pathologic tumor stage | List | Stage of cancer | Stage I, Stage II, Stage III, Stage IIIA |
| days_to_death | List | Number of days the donor patient died after indexing the sample. | 724, 819 |
| subtype | List | Distinct molecular subtypes of tumor. | COC3, COC2, COC1 | 
| tumor_status | List | Description of the sample for tumor status. | TUMOR FREE, WITH TUMOR |
| tumor_stage | List | Depicts the stage of tumor of the donor patient. | Adverse, Intermediate |
| tumor_grade | List | Description of a tumor on the basis of characteristics of tumor tissue and cells under the microscope. | not reported, G2, G3 |
| classification_of_tumor | List | Types of tumor classified on the basis of genomic, proteomic and other molecular analysis. | not reported, Uveal Melanoma, acute myeloid leukemia |
| progression_or_recurrence | List | History of cancer spread or recurrent tumor. | not reported |
| days_to_last_follow_up | List | Number of days after the last follow up with the donor patient. | 1939, 947 |
| history_other_malignancy | List | History of the patient suffering from other cancer type. | [Not Available] |
| history_neoadjuvant_treatment | List | History of patient undergoing neoadjuvant therapy | No, Yes |
| new_tumor_event_dx_indicator | List | New tumor after initial treatment. | YES, No |
| treatment_outcome_first_course | List | Definition of the disease state on the basis of recurrence of tumor. | [Unknown], Complete, Remission/Response |

### 1.5 Molecular Identifiers

Molecular identifiers have been standardized for identifying molecules of interest from various sources. The following are the molecular identifiers used for the data in Omix Atlas.

### 1.5.1 CCLE

| Field    | Type   |
| :---------- | :---------- |
| Mutation | HUGO Gene Symbol |
| Transcriptomics | HUGO Gene Symbol |
| miRNA | miRBase ids |
| Metabolomics | Refmet |
| Proteomics | RPPA Antibody ids |

### 1.5.2 DepMap

| Field    | Type   |
| :---------- | :---------- |
| Mutation | HUGO Gene Symbol |
| Transcriptomics | HUGO Gene Symbol |
| miRNA | miRBase ids |
| Metabolomics | Refmet |
| Proteomics | RPPA Antibody ids |

### 1.5.3 LINCS, GTEX & GEO

| Field    | Type   |
| :---------- | :---------- |
| Transcriptomics | HUGO Gene Symbol |

### 1.5.4 TCGA
| Field    | Type   |
| :---------- | :---------- |
| Mutation | HUGO Gene Symbol |
| Transcriptomics | HUGO Gene Symbol |
| miRNA | miRBase ids |
| Copy Number | HUGO Gene Symbol | 
| Proteomics | RPPA Antibody ids |
| Methylation | CPG island Ids | 

