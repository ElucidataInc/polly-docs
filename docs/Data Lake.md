#Introduction  

With Polly discover you can access public databases which have been curated and stored in the form of data lakes or make data lakes with your own data. These datalakes can be explored and analyzed either through the Polly Notebook interface on Polly or through the Polly Discover application.

Polly Discover consists of the following major components:

*   Data Lake Curation

    Data lakes are reservoirs of information that contain multi-omics data, annotations from publicly available databases, publications, etc. Reservoirs are further segregated into two parts: public data repositories which are curated by Elucidata using public data sources, and private data repositories, where you can add your proprietary data.  
    <br />

*   Data Lake Exploration: 

    To explore a data lake, Polly provides tools such as Shiny applications and Polly Notebooks. These tools enable you to find relevant data by searching for keywords associated with the file name, file metadata, or the contents of a file.  
    <br />

*   On-the-fly Analysis: 

    Once you have narrowed down relevant omics datasets, you can analyze the dataset(s) on the fly using various statistical analyses, displaying intuitive visualizations and allowing you to create a hitlist while analyzing multiple datasets simultaneously.

##Available public data repositories

Public data repositories on Polly consist of processed and curated datasets from various sources. They can be readily used for searching for new datasets or running an analysis on one or more datasets.

*   AML: Microarray and RNA Sequencing datasets for Acute Myeloid Leukemia.
    
*   GBM: Microarray and RNA Sequencing datasets for Gliblastoma Multiforme.
    
*   IBD: Microarray and RNA Sequencing datasets for Inflammatory Bowel Disease.
    
*   GEO: Microarray and RNA Sequencing datasets from Gene Expression Omnbius.
    
*   Single cell Atlas: Single cell RNA Sequencing datasets from Gene Expression Omnibus
    
*   Single cell kidney: Single cell RNA Sequencing datasets from kidney studies.
    
*   GTEX: Normal tissue RNA Sequencing datasets from Genotype-Tissue Expression project
    
*   TCGA: Tumor RNA Sequencing datasets from The Cancer Genome Atlas.
    

Additionally, the public data repositories also consist of publicly available databases that have been curated for annotations. These publicly available databases are currently part of these repositories.

*   HMDB: Pathway information from Human Metabolome Database.
    
*   KEGG: Pathway information from Kyoto Encyclopedia of Genes and Genomes.
    
*   Reactome: Pathway information from Reactome
    
*   GO: Gene Ontology from GO database.
    
*   GWAS: Phenotypic data from Genome-Wide Association Studies Catalogue.
    
*   GENE_INFO: General Gene Information from Human Protein Atlas.

    
#Polly Discover App

##Opening the app

Upon opening the Discover portal on Polly, choose a data repository that you would like to explore. The page should look something like this.

![Polly Discover](../img/Discover/image1.jpg)

After selecting a repository, you’ll be able to see a dashboard and the Polly Discover app icon. The dashboard will give you an overview of the repository, for instance, the number of data sets in the repository, distribution of tissues, organisms, etc. among the datasets.

![Repository Dashboard](../img/Discover/image2.png)

Click on the below icon to start the discover app.

![Discover App](../img/Discover/image3.png)

The app will open and you should see the overview page which contains a brief **description** of the application, it's **scope** and the **usage** as shown below.

![App Description](../img/Discover/image4.png)

##Exploring the data lake

Search for relevant datasets by navigating to the *Dataset Search* tab in the navigation pane to the left. Keyword search can be applied to the following fields:

*   Data Set ID
    
*   Data Set Source
    
*   Description
    
*   Diseases
    
*   Is Public
    
*   Organisms
    
*   Platform
    
*   Tissue
    
*   Year
    

![Search options](../img/Discover/image5.png)

The search will return all datasets that are associated with your search. The result should look like the image below.

![Search results](../img/Discover/image6.png)

The table shown above shows very few columns by default. In order to view the other columns in the table, you can select the fields from **Available Columns** and click on *Show!* button. *Download Selected Dataset* button will let you download the dataset that you have selected on your local system. *Export results to CSV* button will let you download the search result table in the form of a .csv file. Once you have narrowed down the relevant datasets, you can analyze one or more datasets on the fly within the app.

##Analyzing a single dataset

You can analyze a single dataset by selecting the checkbox to the left of the entry in the table. Once you’ve selected the checkbox, click on the *Analyze Data* button below the table description.

![Select a data set](../img/Discover/image7.png)

After clicking the *Analyze Data* button, the app will read the selected dataset and take you to the *Dataset Analysis* tab. Here, you can perform the following analyses:

*   Principal Component Analysis (PCA)
    
*   Pathway Visualization
    
*   Differential Expression
    
*   X2K Analysis
    
*   Enrichr
    
*   GSEA  

![Analyses possible](../img/Discover/image8.png)

##Analyzing multiple datasets

It is also possible to analyze multiple datasets simultaneously. If you select multiple datasets from the search result and click the *Analyze Data* button, you will see a new *Multiple Analyses* Tab, wherein you can select algorithms that you want to run on the selected datasets. While selecting algorithms, a workflow of nodes is formed. These nodes are input blocks for each algorithm for each dataset you selected.

**Total number of nodes = number of datasets selected (n) X number of algorithms selected (m)**

![Multiple analyses](../img/Discover/image9.png)

Upon clicking each node, a pop up with input parameters for each algorithm corresponding to a particular dataset will appear. Adjust these parameters to fit your analysis before executing the algorithms.

![Algorithm input parameters](../img/Discover/image10.png)

Once you have filled all the parameters, you can verify them in the table formed below. Next, execute the algorithms by clicking on the *Execute Algorithms* button.

![Execute algorithms](../img/Discover/image11.png)

After the algorithms have been executed, you’ll be taken to the rules tab, where you can define the rules by which you want to filter the output of the algorithms to get a preliminary hitlist. For instance, after running differential expression, you want to filter by logFC and *p*-value. Along with that you can also add annotations to that hitlist by selecting your source of annotations.

![Rules tab](../img/Discover/image12.png)

Finally, on clicking the *Execute Rules* button, you will be able to get a preliminary hitlist.

![Preliminary hitlist](../img/Discover/image13.png)

#Access through Polly Notebook interface

The Polly Notebook Dockers on Polly have an **internal** python package called ‘discoverPy’ pre-installed, which can be used to search for datasets in the various data repositories.

##Structure of a data repository

A data repository is a collection of different files having different file types. To ensure easy access at a granular level to all datasets a data repository is organized in the following manner. Under this schema, each repository can be considered as a collection of indices which can be used for querying. The discoverPy package can access all indices of a data repository using API endpoints.<!--Structure of a data repository-->

A data repository is a collection of different files having different file types. To ensure easy access at a granular level to all datasets a data repository is organized in the following manner. Under this schema, each repository can be considered as a collection of indices which can be used for querying. The discoverPy package can access all indices of a data repository using API endpoints.  

![Structure of a data repository](../img/Discover/discover_data_repository.png)

Click [here](../Scaling compute/Polly Notebooks) for a detailed documentation about Polly Notebooks.

##Usage

*   Initialize a discover object

This discover object is used to interact with a dataset repository.

<pre><code>from discoverpy import Discover
discover = Discover() 
discover</code></pre>

![Discover object](../img/Discover/image14.png)

*   List all available endpoints for public data repositories along with their index

<pre><code>discover.get_endpoints()</code></pre>

![Endpoints for public data](../img/Discover/image15.png)

*   Add a dataset repository file index

This index contains all file names, their paths in the cloud along with any relevant metadata.

<pre><code>discover.set_dataset_repo("aml_files")</code></pre>

*   Add a sample index

This index contains information at the sample level (columns of a matrix) with any relevant metadata.

<pre><code>discover.set_dataset_repo("aml_gct_metadata")</code></pre>

*   Add a feature index
This index contains information at the feature level (rows of a matrix). Every data value in the matrix is saved in this repository feature wise.

<pre><code>discover.set_feature_repo("aml_gct_data")</code></pre>

After you’ve added the indices for a repository, you can view the discover object

<pre><code>discover</code></pre>

![View the discover object](../img/Discover/image16.png)

Note that the ‘annotation\_repo’ index is added automatically for each repository.

##Querying at the dataset level

To search for datasets, the ‘\_files’ index can be searched using the metadata fields present in it.

*   Get fields present in the index

<pre><code>discover.dataset_repo.get_all_fields()</code></pre>

![Fields present in the index](../img/Discover/image17.png)

*   To get a sense of what values are present in each field, one can view the top n entries. Some generic fields are present for each file.

    *   \_\_bucket\_\_: S3 bucket name
    
    *   \_\_filetype\_\_: Type of file such as pdf, gct etc.
    
    *   \_\_key\_\_: S3 key of the file
    
    *   \_\_location\_\_: Location of the file within data repository

<pre><code>discover.dataset_repo.get_top_n_examples(n = 30)</code></pre>

![Values present in each field](../img/Discover/image18.png)

*   Search for a dataset by keyword in a particular field. Searching for “mll” in the field “description” here.

<pre><code>dataset_query_df = discover.dataset_repo.query_dataset_by_field("description","mll") dataset_query_df</code></pre>

![Search for a dataser](../img/Discover/image19.png)

*   Download file for a dataset by key

<pre><code>discover.download_file_by_key(dataset_query_df["_source.__key__"][14])</code></pre>

##Querying at the sample level

GCT File Format

The datasets in the public repositories are saved as a .gct file. This is a file format in which data can be stored along with the sample metadata. The data values in the actual matrix along with features (genes) are indexed in the ‘\_gct\_data’ index of the repository and the sample metadata is index in the ‘\_gct\_metadata’ of the index of the repository.

![GCT file structure](../img/Discover/gct_file.png)

*   Get fields present in the index

<pre><code>discover.sample_repo.get_all_fields()</code></pre>

![Fields present in the index](../img/Discover/image20.png)

*   To get a sense of what values are present in each field, one can view the top n entries.

<pre><code>discover.sample_repo.get_top_n_examples(n =30)</code></pre>

![Values present in field](../img/Discover/image21.png)

*   Search for samples by keyword in a particular field. Searching for “M1” in the field “fab\_classification\_ch1” here.

<pre><code>fab_df = discover.sample_repo.query_samples_by_field("fab_classification_ch1", "M1", n = 100) fab_df</code></pre>

![Search for samples using keywords](../img/Discover/image22.png)

*   Search for samples by keywords in all fields. This can be used if the field to search for is not known beforehand.

<pre><code>fab_all_fields = discover.sample_repo.query_samples_by_all_fields("M1", n = 100) fab_all_fields</code></pre>

![Search for samples using keywords](../img/Discover/image23.png)

##Querying at the feature level

The matrix of a .gct file contains the actual values for the different features(genes/metabolites). The ‘\_gct\_index’ index of a repository can be queried for features.

*   Get fields present in the index  

<pre><code>discover.feature_repo.get_all_fields()</code></pre>

![Fields present in the index](../img/Discover/image24.png)

*   To get a sense of what values are present in each field, one can view the top n entries.

<pre><code>discover.feature_repo.get_top_n_examples(n=20)</code></pre>

![Values present in field](../img/Discover/image25.png)

The ‘\_\_index\_\_’ column contains the feature name

*   Get values for a particular feature across all samples in all datasets of the repository. Getting values for “NRAS” gene here.

<pre><code>nras_df =discover.feature_repo.get_feature_values("NRAS", n = 1000) nras_df</code></pre>

![Values for a feature](../img/Discover/image26.png)

##Access annotation repositories

The various gene annotation databases can also be accessed through discoverpy. These can be used to get information about a particular gene or a set of genes.

*   Get all annotation databases

<pre><code>discover.annotation_repo.get_annotation_databases()</code></pre>

![Get all annotation databases](../img/Discover/image27.png)

*   Get annotations for a list of genes from a particular database. Getting Reactome pathways for the genes.

<pre><code>discover.annotation_repo.get_feature_annotation('reactome', ['ACTA1', 'AHCTF1', 'AKAP13', 'ATP2C1', 'CDK7'])</code></pre>

![Get annotations for a list of genes](../img/Discover/image28.png)

##Advanced queries

You can also perform more complex queries on multiple fields combining them with boolean logic. Some examples are shown here.

*   Get microarray stem cell datasets which did not involve a knockdown experiment

<pre><code>discover.dataset_repo.query_dataset_by_field_combination(and_fields={"platform":"Microarray", "tissue":"stem cells"}, not_fields={"description":"knockdown"}, n = 50)</code></pre>

![Advanced query example 1](../img/Discover/image29.png)

*   Get samples containing CD34 cells or mononuclear cells do not include de novo samples

<pre><code>discover.sample_repo.query_samples_by_field_combination(or_fields = {"cell_type_ch1":"CD34","cell_type_ch1":"mononuclear"}, not_fields = {"treatment_protocol_ch1":"de novo"}, n = 300)</code></pre>

![Advanced query example 1](../img/Discover/image30.png)

