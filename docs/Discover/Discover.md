# Polly Discover

## Introduction  

With Polly discover you can access public databases which have been curated and stored in the form of data lakes or make data lakes with your own data. These datalakes can be explored and analysed either through the jupyter notebook interface on Polly or through the Polly Discover application.

Polly Discover consists of the following major components:

1.  Data Lake Curation
    
2.  Data Lake Exploration
    
3.  On-the-fly Analysis
    

Data Lake Curation: Data lakes are reservoirs of information that contain multi-omics data, annotations from publicly available databases, publications, etc. Reservoirs are further segregated into two parts: public data repositories which are curated by Elucidata using public data sources, and private data repositories, where you can add your proprietary data.

Data Lake Exploration: To explore a data lake, Polly provides tools such as Shiny applications and Jupyter notebooks. These tools enable you to find relevant data by searching for keywords associated with the file name, file metadata, or the contents of a file.

On-the-fly Analysis: Once you have narrowed down relevant omics datasets, you can analyze the dataset(s) on the fly using various statistical analyses, displaying intuitive visualizations and allowing you to create a hitlist while analyzing multiple datasets simultaneously.

### Available public data repositories  

Public data repositories on Polly consist of processed and curated datasets from various sources. They can be readily used for searching for new datasets or running an analysis on one or more datasets.

1.  AML - Microarray and RNA Sequencing datasets for Acute Myeloid Leukemia.
    
2.  GBM - Microarray and RNA Sequencing datasets for Gliblastoma Multiforme.
    
3.  IBD - Microarray and RNA Sequencing datasets for Inflammatory Bowel Disease.
    
4.  GEO - Microarray and RNA Sequencing datasets from Gene Expression Omnbius.
    
5.  Single cell Atlas - Single cell RNA Sequencing datasets from Gene Expression Omnibus
    
6.  Single cell kidney - Single cell RNA Sequencing datasets from kidney studies.
    
7.  GTEX - Normal tissue RNA Sequencing datasets from Genotype-Tissue Expression project
    
8.  TCGA - Tumor RNA Sequencing datasets from The Cancer Genome Atlas.
    

Additionally, the public data repositories also consist of publicly available databases that have been curated for annotations. These publicly available databases are currently part of these repositories.

1.  HMDB - Pathway information from Human Metabolome Database.
    
2.  KEGG - Pathway information from Kyoto Encyclopedia of Genes and Genomes.
    
3.  Reactome - Pathway information from Reactome
    
4.  GO - Gene Ontology from GO database.
    
5.  GWAS - Phenotypic data from Genome-Wide Association Studies Catalogue.
    
6.  GENE\_INFO - General Gene Information from Human Protein Atlas.

    
## Polly Discover App

### Opening the app

Upon opening the Discover portal on Polly, choose a data repository that you would like to explore. The page should look something like this.

![](./img/image1.png)

After selecting a repository, you’ll be able to see a dashboard and the Polly Discover app icon. The dashboard will give you an overview of the repository, for instance, the number of data sets in the repository, distribution of tissues, organisms, etc. among the datasets.

![](./img/image2.png)

Click on the below icon to start the discover app.

![](./img/image3.png)

The app will open and you should see the overview page shown below. It basically contains a brief **description** of the application, it's **scope** and the **usage**.

![](./img/image4.png)

### Exploring the data lake

Search for relevant datasets by navigating to the *Dataset Search* tab in the navigation pane to the left. Keyword search can be applied to the following fields:

-   Data Set ID
    
-   Data Set Source
    
-   Description
    
-   Diseases
    
-   Is Public
    
-   Organisms
    
-   Platform
    
-   Tissue
    
-   Year
    

![](./img/image5.png)

The search will return all datasets that are associated with your search. The result should look like the image below.

![](./img/image6.png)

The table shown above shows very few columns by default. In order to view the other columns in the table, you can select the fields from **Available Columns** and click on *Show!* button. *Download Selected Dataset* button will let you download the dataset that you have selected on your local system. *Export results to CSV* button will let you download the search result table in the form of a CSV file. Once you have narrowed down the relevant datasets, you can analyze one or more datasets on the fly within the app.

### Analyzing a single dataset

You can analyze a single dataset by selecting the checkbox to the left of the entry in the table. Once you’ve selected the checkbox, click on the *Analyze Data* button below the table description.

![](./img/image7.png)

After clicking the *Analyze Data* button, the app will read the selected dataset and take you to the *Dataset Analysis* tab. Here, you can perform the following analyses:

-   Principal Component Analysis (PCA)
    
-   Pathway Visualization
    
-   Differential Expression
    
-   X2K Analysis
    
-   Enrichr
    
-   GSEA
    

![](./img/image8.png)

### 

Analyzing multiple datasets

It is also possible to analyze multiple datasets simultaneously. If you select multiple datasets from the search result and click the *Analyze Data* button, you will see a new *Multiple Analyses* Tab, wherein you can select algorithms that you want to run on the selected datasets. While selecting algorithms, a workflow of nodes is formed. These nodes are input blocks for each algorithm for each dataset you selected.

**Total number of nodes = number of datasets selected (n) X number of algorithms selected (m)**

![](./img/image9.png)

Upon clicking each node, a pop up with input parameters for each algorithm corresponding to a particular dataset will appear. Adjust these parameters to fit your analysis before executing the algorithms.

![](./img/image10.png)

Once you have filled all the parameters, you can verify them in the table formed below. Next, execute the algorithms by clicking on the *Execute Algorithms* button.

![](./img/image11.png)

After the algorithms have been executed, you’ll be taken to the rules tab, where you can define the rules by which you want to filter the output of the algorithms to get a preliminary hitlist. For instance, after running differential expression, you want to filter by logFC and p-value. Along with that you can also add annotations to that hitlist by selecting your source of annotations.

![](./img/image12.png)

Finally, on clicking the *Execute Rules* button, you will be able to get a preliminary hitlist.

![](./img/image13.png)

## Access through Jupyter notebook interface  

The jupyter notebook Dockers on Polly have an **internal** python package called ‘discoverPy’ pre-installed, which can be used to search for datasets in the various data repositories.

### Structure of a data repository

A data repository is a collection of different files having different file types. To ensure easy access at a granular level to all datasets a data repository is organised in the following manner. Under this schema, each repository can be considered as a collection of indices which can be used for querying. The discoverPy package can access all indices of a data repository using API endpoints.Structure of a data repository

A data repository is a collection of different files having different file types. To ensure easy access at a granular level to all datasets a data repository is organised in the following manner. Under this schema, each repository can be considered as a collection of indices which can be used for querying. The discoverPy package can access all indices of a data repository using API endpoints.  

![](./img/discover_data_repository.png)

### Usage

1\. Initialise a discover object.  
This discover object is used to interact with a dataset repository.

<pre><code>from discoverpy import Discover
discover = Discover() 
discover</code></pre>

![](./img/image14.png)

2\. List all available endpoints for public data repositories along with their index.

<pre><code>discover.get_endpoints()</code></pre>

![](./img/image15.png)


3\. Add a dataset repository file index.  
This index contains all file names, their paths in the cloud along with any relevant metadata.

<pre><code>discover.set_dataset_repo("aml_files")</code></pre>

4\. Add a sample index.

This index contains information at the sample level(columns of a matrix) with any relevant metadata.

<pre><code>discover.set_dataset_repo("aml_gct_metadata")</code></pre>

5\. Add a feature index.  
This index contains information at the feature level(rows of a matrix). Every data value in the matrix is saved in this repository feature wise.

<pre><code>discover.set_feature_repo("aml_gct_data")</code></pre>

After you’ve added the indices for a repository, you can view the discover object

<pre><code>discover</code></pre>

![](./img/image16.png)

Note that the ‘annotation\_repo’ index is added automatically for each repository.

### Querying at the dataset level

To search for datasets, the ‘\_files’ index can be searched using the metadata fields present in it.

1.  Get fields present in the index
    

<pre><code>discover.dataset_repo.get_all_fields()</code></pre>

![](./img/image17.png)


2\. To get a sense of what values are present in each field, one can view the top n entries. Some generic fields are present for each file.

1.  \_\_bucket\_\_: S3 bucket name
    
2.  \_\_filetype\_\_: Type of file such as pdf, gct etc.
    
3.  \_\_key\_\_: S3 key of the file
    
4.  \_\_location\_\_: Location of the file within data repository
    

<pre><code>discover.dataset_repo.get_top_n_examples(n = 30)</code></pre>

![](./img/image18.png)


3\. Search for a dataset by keyword in a particular field. Searching for “mll” in the field “description” here.

<pre><code>dataset_query_df = discover.dataset_repo.query_dataset_by_field("description","mll") dataset_query_df</code></pre>

![](./img/image19.png)


4\. Download file for a dataset by key

<pre><code>discover.download_file_by_key(dataset_query_df["_source.__key__"][14])</code></pre>

### Querying at the sample level

GCT File Format

The datasets in the public repositories are saved as GCT. This is a file format in which data can be stored along with the sample metadata. The data values in the actual matrix along with features(genes) are indexed in the ‘\_gct\_data’ index of the repository and the sample metadata is index in the ‘\_gct\_metadata’ of the index of the repository.

![](./img/gct_file.png)

1.  Get fields present in the index
    

<pre><code>discover.sample_repo.get_all_fields()</code></pre>

![](./img/image20.png)


2\. To get a sense of what values are present in each field, one can view the top n entries.

<pre><code>discover.sample_repo.get_top_n_examples(n =30)</code></pre>

![](./img/image21.png)


3\. Search for samples by keyword in a particular field. Searching for “M1” in the field “fab\_classification\_ch1” here.

<pre><code>fab_df = discover.sample_repo.query_samples_by_field("fab_classification_ch1", "M1", n = 100) fab_df</code></pre>

![](./img/image22.png)

4\. Search for samples by keywords in all fields. This can be used if the field to search for is not known beforehand.

<pre><code>fab_all_fields = discover.sample_repo.query_samples_by_all_fields("M1", n = 100) fab_all_fields</code></pre>

![](./img/image23.png)


### Querying at the feature level

The matrix of a GCT file contains the actual values for the different features(genes/metabolites). The ‘\_gct\_index’ index of a repository can be queried for features.

1.  Get fields present in the index
    

<pre><code>discover.feature_repo.get_all_fields()</code></pre>

![](./img/image24.png)


2\. To get a sense of what values are present in each field, one can view the top n entries.

<pre><code>discover.feature_repo.get_top_n_examples(n=20)</code></pre>

![](./img/image25.png)


The ‘\_\_index\_\_’ column contains the feature name

3\. Get values for a particular feature across all samples in all datasets of the repository. Getting values for “NRAS” gene here.

<pre><code>nras_df =discover.feature_repo.get_feature_values("NRAS", n = 1000) nras_df</code></pre>

![](./img/image26.png)

### Access annotation repositories

The various gene annotation databases can also be accessed through discoverpy. These can be used to get information about a particular gene or a set of genes.

1\. Get all annotation databases.

<pre><code>discover.annotation_repo.get_annotation_databases()</code></pre>

![](./img/image27.png)


2\. Get annotations for a list of genes from a particular database. Getting Reactome pathways for the genes.

<pre><code>discover.annotation_repo.get_feature_annotation('reactome', ['ACTA1', 'AHCTF1', 'AKAP13', 'ATP2C1', 'CDK7'])</code></pre>

![](./img/image28.png)


### Advanced queries:-

One can also perform more complex queries on multiple fields combining them with boolean logic. Some examples are shown here.

1.  Get microarray stem cell datasets which did not involve a knockdown experiment
    

<pre><code>discover.dataset_repo.query_dataset_by_field_combination(and_fields={"platform":"Microarray", "tissue":"stem cells"}, not_fields={"description":"knockdown"}, n = 50)</code></pre>

![](./img/image29.png)


2\. Get samples containing CD34 cells or mononuclear cells do not include de novo samples

<pre><code>discover.sample_repo.query_samples_by_field_combination(or_fields = {"cell_type_ch1":"CD34","cell_type_ch1":"mononuclear"}, not_fields = {"treatment_protocol_ch1":"de novo"}, n = 300)</code></pre>

![](./img/image30.png)

