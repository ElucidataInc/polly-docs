# Release Notes
<!--August 14, 2023-->

<details open>
<summary><font size="+1"><b>August 14, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
 <ul> 
<li>Polly-python:
<ul>
<li>Users can now perform meta-analysis using Polly-python and a shiny app derived from it. 
</li></ul>
</li></ul>
  <br>
      <p class="update-button">Update</p>
      <ul>
<li>Users can now preview reports and PDFs in the workspace GUI.
       <li>Datasets can now be moved instantly for consumption on the UI.
<li>Polly-python:
<ul>
<li>Authentication is getting more standard and secure with bearer token-based authentication.
</li></ul>
</ul>
</details>
<hr>

<!--July 17, 2023-->

<details open>
<summary><font size="+1"><b>July 17, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
 <ul> 
    <li>Cell-type distribution numbers are available for single-cell datasets on the revamped details page under the Cell Type Visualization tab. The Metadata Table tab contains a schema-driven and customizable metadata table that can be used to interact with cell-level metadata, while the Dataset Overview tab provides general dataset details.
    <li>The search results on OmixAtlases have been enhanced, and the keywords found will now be ranked higher. This is because the search bar is equipped with a new NER model that can recognize disease, tissue, drugs, and cell lines.
<li>GEO accession numbers and PubMed IDs will be highlighted in yellow in the search results when searching for datasets using them. Other text keywords will be highlighted with bold text as they appear in the different description tabs on the cards.
       </ul>
  <br>
      <p class="update-button">Update</p>
      <ul>
<li>Switching between card view and table view will be faster than before.
<li>Polly-python:
<ul>
<li>The `query_metadata` function has been optimized for memory utilization.
<li>The `identify_cohorts` function has been improved. It will:
<ul>
        <li>Give distribution of factors to users
<li>Show number of samples in all cohorts
<li>Show users how to plot sunburst with custom columns as code in docs
<li> Print a message when user calls the function with 1 sample and the sunburst is empty
        </ul>
<li>polly-python installation will be smoother than before as it has been tested for compatibility for python versions >3.6.
</li></ul>
</ul>
</details>
<hr>
<!--June 19, 2023-->

<details open>
<summary><font size="+1"><b>June 19, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>A product tour is available in the help section of the navigation bar on the left side. A series of popup messages highlighting the different features and steps will be offered to the users for their onboarding on Polly.
    <li>Users now have the choice to log out of all currently open sessions when they change their passwords. 
       </ul>
       <br>
      <p class="update-button">Update</p>
  <ul>
         <li>Users will now be able to access the datasets they had shortlisted in the Public Data Atlases even after logging out because the datasets in the shortlist will still be available after the user logs out of the current session.
         <li>The total number of datasets in the Public Data Atlases that are available is now reflected accurately in the dataset number on the OmixAtlas homepage.
         <li>When applications are opened, more precise timing and the dataset size are displayed.
       <li>The quality of data with respect to differentiating empty fields has been improved by ensuring consistency between ‘Normal’ inputs in the disease field and 'none' inputs in other fields. 
       </li>
       </ul>
</details>
<hr>
<!--May 22, 2023-->

<details open>
<summary><font size="+1"><b>May 22, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>The Details page of datasets in the Bulk RNAseq OmixAtlas has been revamped to offer a cleaner interface with information divided into three sections:
           <ul>
                  <li>Dataset Overview: This section provides basic details about the dataset such as PubMed ID, link to the source, authors, summary, etc. 
                  <li>Metadata Table: A table with cleaner and more harmonized metadata columns compared to the data source.
                  <li>Metadata Charts: A section where you can create interactive sunburst plots with up to 4 metadata fields plotted at a time. By default, experimental factors are represented in the sunburst plot.
           </ul> 
    <li>Experimental variables that vary within the samples of a dataset are available as a separate field on Polly called `experimental factors`. This field will contain a list of variables. The list of variables can be found on the table view of the Bulk RNAseq OmixAtlas, in the metadata tables of datasets, and in the dataset overview section of the Details page of datasets. 
    <li>Interactive filters and more relevant graphs are available in the summary pages of OmixAtlases.
    <li>The metadata fields available on Polly’s Bulk RNAseq OmixAtlas have been given cleaner, more intuitive, and more harmonized names to make them easier to understand. 
           <ul>
                  <li>All the source metadata at the sample level will be available under cleaner and harmonized column names. Eg. `characteristics_ch1_2_treatment` → `treatment` , `characteristics_ch1_1_genotype` → `genotype` , `strain_ch1` becomes `strain`
                  <li>Merged columns where information was spread across two or more columns. The merged information is present in one unique column name. Eg. `characteristics_ch1_0_cell_type` , `characterstics_ch1_1_cell_type` → `cell_type`
           </ul>
       </ul>
       <br>
      <p class="update-button">Update</p>
  <ul>
         <li>When users return to the Monitoring Dashboard, the filters that they had previously set will be retained.
         <li>When users attempted to transfer datasets from the public data atlas to the user's data atlas, the ingestion process used to be triggered 15 minutes after the attempt. This gap has been removed, and the ingestion process begins immediately after the transfer has been triggered.
         <li>OmixAtlas landing pages will load faster than before, almost instantly.
         <li>The majority of the single cell datasets load on Cellxgene within 25-40 seconds, much faster than before.
       </ul>
</details>
<hr>
<!--April 24, 2023-->

<details open>
<summary><font size="+1"><b>April 24, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>With cell-type ontology recommendations on the OmixAtlas filters, users receive matching and related cells for the cell types they enter. These recommendations contain cells that are hypernyms, hyponyms, or a part of the given cell type. If the user enters the name of a tissue, all the cells that are part of that tissue will be recommended to the user.
    <li>New datasets on GEO will be added to the Bulk RNASeq OmixAtlas every two weeks as they come on to GEO.
    <li>The Monitoring Dashboard has been enhanced with optimized filters, and more accurate statuses for runs.
    <li>While querying on polly-python, users can now query based on ‘views’ available in the OmixAtlas. This will help generate query responses specific to the dataset source and/or data type mentioned in the schema.
       </ul>
       <br>
      <p class="update-button">Update</p>
  <ul>
         <li>The failure issues related to transferring datasets from the Public Data Atlases to the User Data Atlases have been fixed.
       </ul>
</details>
<hr>
<!--March 27, 2023-->

<details open>
<summary><font size="+1"><b>March 27, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>The Bulk RNA-Seq and Single Cell OmixAtlases have a shortlist page that can be accessed using the shortlist icon on the top left corner of the View Datasets screen. Datasets that users want to shortlist before buying will appear here.
    <li>Multiple datasets can be selected for shortlisting in one go using the check box feature on each dataset both in Card View and Table View of OmixAtlases. 
    <li>Academic users can sign up directly on Polly, for the Play with Polly project, using their .edu or .ac email addresses.
    <li>Polly-python supports data ingestion of .CSV files which may include: datasets, metadata files, analyses results and so on.
    <li>The metadata fields ‘year’ and ‘authors’ have been added to all the datasets in the Bulk RNAseq OmixAtlas.
    <li>Polly is now SOC-2 certified. This means that users' data is secure and protected against data breaches. This certification also overlaps with other industry standards such as ISO 27001, HIPAA, etc.
       </ul>
       <br>
      <p class="update-button">Update</p>
  <ul>
         <li>UI improvements have been made to the side navigation bar on Polly, the OmixAtlas page and the OmixAtlas cards.
         <li>Users will be notified if they are using an older version of polly-python. 
       </ul>
</details>
<hr>
<!--February 27, 2023-->

<details open>
<summary><font size="+1"><b>February 27, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>When users search for keywords on the OmixAtlas search bar, they will also get results that match all the available sample-level metadata.   
    <li>Abstracts are now available for datasets in the the Bulk RNA-seq OmixAtlas. 
    <li>.gct files can now be directly downloaded from Phantasus when accessed through the user’s Destination OmixAtlas. 
    <li>Restricted files, that can only be accessed by users who have access to the workspace in which the file is present, can now be downloaded.
       </ul>
       <br>
       <p class="Polly-python-button">Polly-python</p>
  <ul>
    <li>The move_data function can be used to swiftly move data from the source to the destination OmixAtlases without having to download and re-upload the data.  
    <li>The entire sample level metadata for a dataset can be obtained using the get_metadata function. 
    <li>To download datasets easily, users can use the new download_dataset() function.
    <li>Upon changing the schema of an OmixAtlas, users will be alerted of the impact of the specific change. 
                </ul>
<br>
       <p class="update-button">Update</p>
  <ul>
         <li>On **polly-python**, upon changing the schema of an OmixAtlas, users will be alerted of the impact of the specific change.
         <li>On the OmixAtlases, the performance of omixatlas summaries, dataset details, datasets search and search filters has been improved to facilitate a smoother experience. 
       </ul>
</details>
<hr>
<!--January 30, 2023-->

<details open>
<summary><font size="+1"><b>January 30, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>For base dockers, python, R and pollyglot, memory & CPU utilization will be displayed within the notebooks.  
    <li>Users will be informed when they utilize 70% of the memory in their notebooks.
    <li>The notebook loading page will show the time it will take to load notebooks 
    <li>Drug ontology recommendation is now available on the OmixAtlas UI. Structurally similar drugs will be recommended when users search for drugs and on searching for genes, drugs that target the gene will be suggested.
    <li>OmixAtlas cards on the OmixAtlas home page have tags that differentiate source and destination OmixAtlases.  
    <li>The Single-Cell OmixAtlas has been created with 1840 cell-type curated datasets from multiple sources, including high value publications.
    <li>Manually curated labels for age, gender, donor and sampling site have been added to the Single-Cell OmixAtlas.
    <li>Files, notebooks and analyses can be pinned on workspaces for easy accessibility. 
    <li>The link_report_url function on polly-python can be used to link a dataset_id to a URL. 
                </ul>
<br>
       <p class="update-button">Update</p>
  <ul>
         <li>Reports linked to datasets can now be previewed instead of being downloaded and then viewed.
         <li>The schema is automatically validated before inserting a new schema or updating the existing schema of an OmixAtlas.
         <li>In the download metadata function, users now have the option to select the keys to be field_name or original_name. 
       </ul>
</details>
<hr>
<!--January 6, 2023-->

<details open>
<summary><font size="+1"><b>January 6, 2023</b></font></summary>
<br>
       <p class="Releases-button">Releases</p>
  <ul> 
    <li>For HTML files in the workspaces, users can download the dependencies (images, folders, CSS files, linked HTML files, etc.) via the UI. 
    <li>Deleted files/folders from workspaces can be retrieved for up to 15 days.
    <li>Drug ontology recommendation will be available via polly-python. Structurally similar drugs will be recommended when users search for drugs and on searching for genes, drugs that target the gene will be suggested. 
    <li>The following fields have been updated for all GEO single-cell datasets:
           <ul>
                  <li>Datasets level: Tissue, Disease, Cell Line, Organism, Gene, Strain, Drug.
                  <li>Sample Level: Tissue, Strain, Disease, Drug, Genotype, Cell Line, Control Perturbation.
           </ul>
    <li>Sample-level cancer stage and dataset-level treatment and clinical labels were added to the new BulkRNASeq OmixAtlas with >42k datasets (Kallisto raw counts). 
    <li>Users now have an option to move datasets from the Bulk RNA-Seq and Single Cell OmixAtlases to their organization’s OmixAtlases by themselves.
    <li>Users can request for additional services such as curation and QC directly from the OmixAtlas UI.
                </ul>
<br>
       <p class="update-button">Update</p>
  <ul>
         <li>The workspace search is faster, offering better results based on the user’s query.
         <li>New single-cell datasets from Broad SC Portal, ExpressionAtlas and Covid19 Atlas have been added to the Single-Cell OmixAtlas.
         <li>The speed of dataset ingestion for Bulk RNA-seq data has been improved 30x. 100 datasets are currently being added per day.
         <li>Higher anndata version (>=0.8) is supported while indexing .h5ad files.
         <li>Arrays can now be stored in sample level metadata (.gct files).
         <li>The visibility of sample level metadata in the Details Page has been improved for understanding a dataset better.
       </ul>
</details>
<hr>
<!--December 5, 2022-->

<details open>
<summary><font size="+1"><b>December 5, 2022</b></font></summary>
<br>
       <p class="releases-button">Releases</p>
  <ul> 
    <li>.h5seurat files are now supported on OmixAtlases. 
    <li>3 new fields have been curated for 101k GEO datasets on the platform - cancer stage, chemical treatment and clinical labels.
    <li>On the OmixAtlas, the filter result counts for each category are visible and clickable to enable easy navigation to relevant datasets. For instance, on filtering diseases, users will get 18414 normal samples, 2063 neoplasms and 1374 obesity samples which users can access by just clicking on the respective filter entries.
    <li>The OmixAtlas table view can be expanded now since the metadata fields visible on the view can be customized through the schema. This offers flexibility for locating useful datasets. 
                </ul>
<br>
       <p class="update-button">Update</p>
  <ul>
         <li>Ontology recommendation on OmixAtlases has been improved. Now, the recommendations will be more relevant, accurate and with more unique suggestions. 
         <li>40k RNA-seq datasets from GEO were added to a new OmixAtlas.
       </ul>
</details>
<hr>
<!--November 7, 2022-->

<details open>
<summary><font size="+1"><b>November 7, 2022</b></font></summary>
<br>
       <p class="new-button">New</p>
  <ul>
    <li>Using the cost dashboard on the Polly Admin panel, an organization admin can track the compute cost for every user based on several parameters - machine type, usage of apps, jobs and notebooks, etc.
    <li>Users can now reorganize the notebooks within workspaces, at their convenience, by creating and cloning notebooks inside folders.
    <li>Users will receive notifications detailing the data sources and types of the new datasets ingested into the OmixAtlas that they have subscribed to.
    <li>The free text search bar on OmixAtlases enables search for keywords across dataset-level metadata fields. E.g., title, description, summary, tissue, drug, disease, etc., depending on the schema of the OmixAtlas.
In this release, it supports even more precise and advanced search operations using logical operators:      
      <ul>
        <li>And &  
        <li>Or |
        <li>Not ~
        <li>Group ()
        <li>Exact match ""
                </ul>
            </ul>
<br>
       <p class="update-button">Update</p>
  <ul>
         <li>Using the schema management module of Polly-python, users can customize the following on the OmixAtlas: 
                <ul>
                       <li>columns for table view 
                       <li>filters
                       <li>availability of ontology recommendations 
                       <li>search fields 
                </ul>
                <li>Curated metadata from the curation app can also be exported as a .
                       file apart from .json.  
                <li>New users can accept the End User License Agreement (EULA) directly on the Polly platform during their first login. 
       </ul>
</details>
<hr>
<!--October 28, 2022-->

<details open>
<summary><font size="+1"><b>October 28, 2022</b></font></summary>
<br>
       <p class="new-button">New</p>
  <ul>
    <li>Users can now run jobs, cancel jobs and fetch status of jobs.
    <li>There is a new function for schema validation.
    <li>Users can download dataset level metadata.
    <li>Users can geenrate a merged GCT file from a cohort.
  </ul>
</details>
<hr>
<!--October 10, 2022-->

<details open>
<summary><font size="+1"><b>October 10, 2022</b></font></summary>
<br>
       <p class="new-button">New</p>
  <ul>
    <li>A card view layout page is available where you can sort and filter workspace cards. Tags can be attached with workspaces for better findability.
    <li>Users can mark workspaces as favorites and also watch workspaces to start receiving notifications specific to them. They can be sorted by their creation date for easier access.
    <li>The curation app can be accessed through the side navigation bar on Polly. Users can now switch between reviewer view and curator view directly from the UI.
    <li>With the newly deployed ingestion dashboard, users can track all ingestion runs made to OmixAtlas in real-time, view logs in the event of a failure, and list the completed ingestion runs.     
      <ul>
        <li>With the newly deployed ingestion dashboard, users can track all ingestion runs made to OmixAtlas in real-time, view logs in the event of a failure, and list the completed ingestion runs. 
        <li>Track running time and logs is possible for failed/rescheduled jobs too.
                </ul>
    <li>Schedule compute jobs and apps(notebooks and shiny) on separate nodes to optimize/reduce cost for both groups and increase stability. 
    <li>Support for multiple docker images in the same Nextflow workflow is available.
    <li>The drug labels for Tier 2 datasets on GEO OA are now more accurate and relevant. 
           <ul>
                  <li>All datasets and samples are now annotated with Pubchem identifiers (a change from CHEBI) and encompass a wider range of drug classes (including monoclonal antibodies).
                  <li>The more accurate disease labels for Tier 2 datasets on GEO OA will lead to an improved cohorting and search experience.
           </ul>
            </ul>
         <b>Polly-python:</b>
       <ul>
          <li>The curation library is integrated with polly-python which will enable users to recognize entities in a given text, standardize them and tag them in a text with standardized nomenclature/ontology.
          <li>Schema management-related functions are now upgraded to enable users to update and replace schema. Update schema is to be used to make minor edits in the existing schema, and replace schema is to be used for completely replacing the existing schema of an OmixAtlas. 
          <li>Schema validation functions have been released to enable users to validate the schema prior to inserting/updating or replacing the schema.
          <li>Results of query_metadata are now sorted alphabetically across the columns to improve the UX of browsing in the data frame
          <li>The following file types/extensions can be ingested - biom, zip, tar.gz, gct.bz, vcf.bgz, fcs, fs.
       </ul>
     <br>  
       <p class="update-button">Update</p>
  <ul>
         <li>Data Addition: 
                <ul>
                       <li>30k datasets were added from HugeAMP and OpenGWAS. 
                       <li>HTAN single cell studies were added.
                </ul>
       </ul>
</details>
<hr>
<!--September 19, 2022-->

<details open>
<summary><font size="+1"><b>September 19, 2022</b></font></summary>
<br>
    <p class="update-button">Update</p>
  <ul>
    <li>.nf files can now be uploaded to workspaces.
    <li>New machines are available to run jobs:-
      <ul>
        <li>mix5xlarge: 16 vCPUs, 512.0 GiB RAM  & 2TB storage
        <li>mix6xlarge: 24 vCPUs, 768.0 GiB & 2TB storage
        <li>mix7xlarge: 64 vCPUs, 1024.0 GiB & 3TB storage
                </ul>
    <li>A new version of Polly-CLI is available. It includes performance improvements and bug fixes related to job logs. 
    <li>If the OmixAtlas has multiple sources and/or datatypes, user can now put a schema specific to all the sources and/or datatypes. They will also be able to query the table as specific to source and/or datatype as defined in the schema.
    <li>1500 curated CPTAC datasets have been added to Polly.
    <li>Cohort creation is now enabled for all types of Public and Enterprise OmixAtlas
    <li>Cell line ontology recommendations are now available on polly-pyhton and the front-end of the OmixAtlases on Polly.
    <li>Card View on the OmixAtlas now allows sorting via recency of data ingestion.
    <li>A View Details Page is now available for the new Expression Atlas datasets.
      <li><b>Polly-python:</b>
       <ul>
          <li>The data matrix within h5ad files (single cell data) can now be queried.
          <li>The curation library is integrated with polly-python which will enable the users to recognise entities in a given text, standardize entity and tag the entities in a text with standardized nomenclature/ontology.
          <li>The replace_schema feature allows users to replace the existing schema entirely. 
          <li>Schema validation functions have been released to enable users to validate the schema prior to inserting/updating or replacing the schema.
          <li>Older reports linked to datasets in any OmixAtlas can now be deleted and newer ones can be linked to the same datasets.
          <li>The rows and columns that appear in the results of query_metadata are now sorted alphabetically to make browsing through the dataframe easier.
          <li>Ingestion of the following file types/extensions is enabled in OmixAtlases - biom, zip, tar.gz, gct.bz, vcf.bgz, fcs, fs.
  </ul>
</details>
    
<hr>

<!--September 5, 2022-->

<details open>
<summary><font size="+1"><b>September 5, 2022</b></font></summary>
<br>
   <p class="new-button">New</p>
  <ul>
    <li>Microbiome data is now supported on OmixAtlases. The data can be searched and queried on Polly and can be downloaded in the form of a BIOM file.
    <li>‘Cancer Stage’ is the newest metadata field and it has been curated for  974 datasets on GEO.
    <li>Reports present in workspaces can now be linked to datasets in OmixAtlases. 
                </ul>
  
  <p class="update-button">Update</p>

<ul>
<li>‘Select All’ is now an option for the filter results while searching and filtering datasets on the OmixAtlas UI. 
<li>Identify datasets that do not contain data matrices on the OmixAtlases through the new ‘Metadata-only’ labels.
<li>Send requests to obtain data matrices for specific datasets. 
  </ul>
</details>
    
<hr>

<!--August, 2022-->

<details open>
<summary><font size="+1"><b>August, 2022</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>The Nextflow integration with Polly CLI allows users to run any Nextflow bioinformatics pipeline with parallel processing on Polly’s scalable computational infrastructure.
    <li>OmixAtlases are now available as shareable links through Workspaces.
    <li>797 clinical datasets from MIMIC are now available on Polly, curated for four metadata fields “Drug, Dose, Frequency, and Strength.” 
    <li>With the recent Polly Python version release, users can now:
            <ul>
              <li>Copy files/folders in workspaces from one workspace to another. 
              <li>Add datasets to and delete datasets from OmixAtlases.
              <li>Link and fetch reports to a dataset on any OmixAtlas. 
              <li>Get auto-generated metadata summaries for datasets present on the GEO OmixAtlas by giving the GEO accession ID as an input. This helps to improve findability and estimate the relevance of the dataset.
              <li>Cell line recommendations are now available, to select multiple related cell lines. Users can start the search with a disease, tissue or cell line and receive recommendations for related or matching cell lines. 
              <li>The recommend functionality is available for disease, tissue and cell lines at sample level metadata queries as well. 
     </ul>
              <li>New datasets from the following repositories have been added to Polly OmixAtlases in the last month:
            <ul>
          <li>DepMap OmixAtlas  -  4312
          <li>GEO OmixAtlas  -  922
          <li>Single-Cell OmixAtlas  -  240
</ul>
                </ul>
  
  <p class="update-button">Update</p>

These are the major Polly Python updates to existing functionalities: 
<ul>
<li>The complete schema for tables in an OmixAtlas can be fetched in the form of a dictionary.
<li>Schema for feature level metadata is now available.
<li>Schema for single cell and GWAS data can be retrieved.
<li>Ontology recommendations for sample-level queries have been enabled.
  </ul>
</details>
    
<hr>
  

<!--March, 2022-->

<details open>
<summary><font size="+1"><b>March, 2022</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Users can now host docker based applications and add their own notebook environments using Polly CLI.</li>
    <li>Users can launch opensource notebooks from GitHub directly on Polly’s compute environment.</li>
    <li>Users can save dataset from OmixAtlas to workspace as well as upload/download files and folders to/from workspace using polly-python.</li>
    <li>Users can filter the schema in an Omixatlas specific to the source & data_type using polly-python.</li>
    <li>Users can access the installed version of polly-python inside python shell or jupyter notebook cell.</li>
    <li>Users can create cohorts for TCGA Transcriptomics and Mutation data using polly-python.</li>
    <li>Users can now create an Omixatlas using polly-python.</li>
    <li>Users can query datasets, samples and features on polly-python across multiple OmixAtlases at once. Find examples here.</li>
    <li>Recommended disease ontologies are displayed when a user queries disease field in OmixAtlas.</li>
    <li>Users can now copy/move files (except analysis files) or folder across workspaces and folders on Polly frontend.</li>
    <li>Users can now launch notebooks situated within folders or sub-folders.</li>
    <li>Users can also filter folders from workspace contents.</li>
  </ul>
  <p class="update-button">Update</p> 
  <ul>
    <li>polly-python users can now access schema functions via both repo_id and repo_name.</li>
    <li>polly-python users can easily convert .gct file format to .maf file format in TCGA and cBioportal repositories using a file format converter function.</li>
    <li>211 Single Cell datasets were added to OmixAtlas with cluster-level cell type annotation.</li>
    <li>For 3.5k GEO datasets platform field was updated on Polly.</li>
    <li>Users can now also preview tsv files along with other files types on Polly.</li>
    <li>Cell types were added for fetal single cell atlas in Polly.</li>
  </ul>
</details>
    
<hr>

<!--January, 2022-->

<details open>
<summary><font size="+1"><b>January, 2022</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Users can create workspaces and fetch list of workspaces using polly-python.</li>
    <li>Change in authentication process - Until now, the users had to authenticate each class separately. In this release, a global authentication mechanism has been developed using which users can authenticate multiple classes (such as OmixAtlas, Workspaces) using a single authentication step.</li>
    <li>Users can publish personal notebook on Polly workspace environment and analyse vcf files using Hail docker.</li>
    <li>Preview of all standard file types is available on Polly UI. Different file types like xls/xlsx, pdf, html, csv, png/jpg/gif, ipynb can be opened directly on Polly without a third party service.</li>
  </ul>
  <p class="update-button">Update</p> 
  <ul>
    <li>2792 datasets, 74713 samples have been annotated for tissue and cell line tags.</li>
    <li>Strain has been added as a queryable field at dataset-level for GEO datasets.</li>
    <li>A bug in Sort by Relevance when searching over description in OmixAtlas table view was fixed. Other sort related bugs have also been fixed.</li>
    <li>360k datasets from the UK Biobank were added on Polly.</li>
    <li>149k Immport lab datasets were added on Polly.</li>
    <li>59k RCSB datasets were added on Polly.</li>
    <li>Cell type curated for 229 Single Cell datasets were added to OmixAtlas.</li>
    <li>Sample level age labels were added to all datasets for TCGA and GEO</li>
  </ul>
</details>
    
<hr>

<!--November, 2021-->

<details open>
<summary><font size="+1"><b>November, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Workspace contents like files, analyses, reports can now be sorted in workspaces based on name, created date, last modified, author, and type.</li>
    <li>15,000 Microarray and WES (Whole Exome Sequencing) datasets from PPMI have been curated and added to Polly.</li>
    <li>Dataset Overview (containing Title, Publication, Abstract, Tags for the Data, sample information as summary plots and table, processed data as a table) for every dataset can be viewed using the “View Details” option beside datasets in OmixAtlases.</li>
    <li>The Datasets (gct, h5ad, vcf) files can be downloaded from the Options Menu in the Card view or from the View Details Page.</li>
    <li>Request for Dataset option is available at multiple places within the OmixAtlases.</li>
  </ul>
  <p class="update-button">Update</p> 
  <ul>
    <li>Resolved a bug in El-MAVEN instance termination.</li>
    <li>Curation information of 45,000 datasets from gnomAD and DepMap has been updated.</li>
  </ul>
</details>
    
<hr>

<!--October, 2021-->

<details open>
<summary><font size="+1"><b>October, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>In addition to Liver OmixAtlas, polly now has GDC, GEO, cBioPortal, PharmacoDB, LINCS and Metabolomics OmixAtlases.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added 17500 datasets on Polly</li>
    <li>Genotype, age and gender annotation were added to 3,900,000 GEO samples</li>
  </ul>
</details>

<hr>

<!--September, 2021-->

<details open>
<summary><font size="+1"><b>September, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>New compute machines Mi5xlarge (32 vCPU, 250GB RAM), Mi6xlarge (64 vCPU, 500GB RAM), Mi7xlarge (64 vCPU, 970GB RAM), GPUsmall (1 GPU, 8 vCPU, 60GB RAM) and GPUmedium (4 GPU, 32 vCPU, 240GB RAM) were added to EL-MAVEN.</li>
    <li>Introduced Polly Files (beta version), a desktop application for transferring files between computer and Polly Workspaces.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Resolved a bug that caused an error in SQL query on “’” expansion.</li>
    <li>gnomAD was enriched with 96,000 new datasets of WES and WGS type.</li>
    <li>Fixed table view bugs and enhanced UI on OmixAtlas.</li>
    <li>New datasets totalling 48,000 were added to Immport, HPA, CPTAC and GTex.</li>
    <li>Auto-curated tags totalling 770,000 were added to polly datasets.</li>
    <li>A bug affecting folder deletion was fixed.</li>
  </ul>
</details>

<hr>

<!--August, 2021-->

<details open>
<summary><font size="+1"><b>August, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced Polly Python Library facilitating powerful search capabilities across dataset, sample, and feature level metadata on any computational environment through code.</li>
    <li>Introduced “View Only Access” on Polly Workspaces – an enterprise grade permission giving more control to admins. </li>
    <li>Enabled Voila Dashboards within Polly Notebooks.</li> 
    <li>Introduced an application resource monitor on EL-Maven, enabling users to monitor the progress of a job and make decisions about requirement of a bigger machine.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Over 155,000 datasets were added to LINCS OmixAtlas on Polly.</li>
  </ul>
</details>

<hr>


<!--July, 2021-->

<details open>
<summary><font size="+1"><b>July, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Liver OmixAltas released.</li>
  </ul>
</details>

<hr>


<!--June 18th, 2021-->

<details open>
<summary><font size="+1"><b>June 18th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Added 12,200 new curated transcriptomics and single cell datasets to various Data Lake.</li>
  </ul>
</details>

<hr>

<!--June 4th, 2021-->

<details open>
<summary><font size="+1"><b>June 4th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Made Dual Mode Visualizatopn application and Untargeted Pipeline Application lighter for heavy datasets to avoid memory leakage.</li>
    <li>Added 3.350 new curated transcriptomics and single cell datasets to various Data Lake.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Resolved issue in Lipidomics, Dual Mode and Polly El-MAVEN Applications.</li>
  </ul>
</details>

<hr>


<!--May 21st, 2021-->

<details>
<summary><font size="+1"><b>May 21st, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced Google Slide Intergration with Polly Notebooks.</li>
     <li>Added Reporting feature using Markdown in Dual Mode Visulaization Application.</li>
    <li>Added 23.350 new curated transcriptomics and single cell datasets to various Data Lake.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added Scree Plot under Quality Check in Dual Mode Visualization Application.
    <li>Resolved forgot password issue.</li>
  </ul>
</details>

<hr>



<!--May 7th, 2021-->

<details>
<summary><font size="+1"><b>May 7th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>El-MAVEN latest beta version is now available on Polly.</li>
    <li>Added 25.120 new curated transcriptomics and single cell datasets to various Data Lake.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added two-way ANOVA capability in Dual Mode visualization Application along with combining multiple conditions or cohorts while performing differential expression.</li>
  </ul>
</details>

<hr>

<!--April 23rd, 2021-->

<details open>
<summary><font size="+1"><b>April 23rd, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced an Admin Dashboard to provide account administrators the convenience to manage their accounts.</li>
    <li>Added 10,120 new curated transcriptomics and single cell datasets to various Data Lake.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>HTML as a data file can now be opened through Workspace itself.</li>
    <li>Resolved an issue to make Dropbox work seamlessly with Workspaces.</li>
  </ul>
</details>

<hr>

<!--April 9th, 2021-->

<details>
<summary><font size="+1"><b>April 9th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Created Shiny and Studio applicationS for feature level search of GEO Datasets.</li>
    <li>Added 11,177 new curated transcriptomics datasets to TCGA Data Lake.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>PDF as a data file can now be opened through Workspace itself.</li>
    <li>Resolved an issue to make Google Drive work seamlessly with Workspaces.</li>
  </ul>
</details>

<hr>


<!--March 26th, 2021-->

<details>
<summary><font size="+1"><b>March 26th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>TEDDY (The Environmental Determinants of Diabetes in the Young) and DEPMAP (Dependency Map) Data Lakes have been added on Polly.</li>
    <li>Introduced option to directly export data to the workspace from a Studio Preset.</li>
    <li>Added 37,177 new curated datasets corresponding to various omics to different Data Lakes.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated Polly Login User Interface.</li>
    <li>Added additional filters to TEDDY Data Lake.</li>
    <li>Resolved issue with app hosting infrastructure to increase stability of apps for better user experience.</li>
  </ul>
</details>

<hr>






<!--March 12th, 2021-->

<details>
<summary><font size="+1"><b>March 12th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced Docker building feature on <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html#docker-management">Polly CLI</a> which enable users to build dockers, check their build status and logs and push dockers to Polly.</li>
    <li>Added 11,470 new curated transcriptomics and single cell datasets to different Data Lakes.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Better accessibility to datasets within OmixWiki with accessibility to metadata filtering options.</li>
  </ul>
</details>

<hr>


<!--February 26th, 2021-->

<details>
<summary><font size="+1"><b>February 26th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced the functionality that enables the users to host their own application on Polly by using <a href="https://docs.elucidata.io/Apps/Host%20Apps.html">Polly CLI</a>.</li>
    <li>Enabled feature level querying for GEO Data Lake.</li>
    <li>Added Genomics docker for variant calling and annotation.</li>
    <li>Added a new notebook environment for Genomics Variant Analysis.</li>
    <li>Enabled partial string search for dataset id in the search bar.</li>
    <li>Added 20,096 new curated transcriptomics and single cell datasets to different Data Lakes.</li>
  </ul>
</details>

<hr>



<!--February 12th, 2021-->

<details>
<summary><font size="+1"><b>February 12th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced the <a href="https://status.polly.elucidata.io/">status page</a> for real time updates on Polly’s status, downtime, incidents, and maintenance.</li>
    <li>Added auto-run feature for selected Studio Presets.</li>
    <li>Enabled component updating and versioning by component creator.</li>
    <li>Added 11,580 new curated transcriptomics datasets to GEO and LINCS Data Lakes.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated the UI of visualization dashboard of <a href="https://docs.elucidata.io/Apps/Data%20Studio/Data%20Studio.html">Data Studio</a> for better visibility.</li>
    <li>Updated all notebook dockers with the latest version of discoverpy (0.0.10).</li>
    <li>Added finer error and warning messages to <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">CLI</a>.</li>
    <li>Removed the 1000 row limit on query results in <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">CLI</a>.</li>
  </ul>
</details>

<hr>


<!--January 29th, 2021-->

<details>
<summary><font size="+1"><b>January 29th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Public sharing of the reports created within any Studio session is now available on Polly.</li>
    <li>Added 14,727 new curated transcriptomics and metabolomics datasets with 9,513 transcriptomics datasets being added to the LINCS Data Lake.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added specific error message to indicate presence of multiple groups with the same compound name in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow</a>.</li>
    <li>Added specific error message in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow</a> if isotopologues of the compound are spread over different metagroups in El-MAVEN output.</li>
  </ul>
</details>

<hr>


<!--January 15th, 2021-->

<details>
<summary><font size="+1"><b>January 15th, 2021</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>GTEx Correlation and Enrichment Analysis preset is now available which can be used to identify enriched pathways based on the gene correlations.</li>
    <li>Added TraceFinder Downstream Analysis preset with additional feature of translating the analytical insights into shareable dashboards.</li>
    <li>Added 1,836 new curated transcriptomics and proteomics datasets to different Data Lakes.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Enabled use of retention time information for metabolite identification and updated <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Untargeted%20Pipeline.html">Untargeted Pipeline</a> library to handle already identified metabolities.</li>
  </ul>
</details>

<hr>


<!--January 1st, 2021-->

<details>
<summary><font size="+1"><b>January 1st, 2021</b></font></summary>
<br>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Untargeted%20Pipeline.html">Untargeted Pipeline</a> to be compatible with El-MAVEN's peakML output.</li>
  </ul>
</details>

<hr>

<!--December 18th, 2020-->

<details>
<summary><font size="+1"><b>December 18th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>LINCS(Library of Integrated Network-Based Cellular Signatures) repository with 19,520 curated datasets has been added in <a href="https://docs.elucidata.io/Data%20Lake.html">Data Lake</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added ANOVA Test and updated Limma Test with extra filters for volcano plot and Heatmap for the differentially expressed results in the <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html">Dual Mode Data Visulaization</a>.</li>
  </ul>
</details>

<hr>

<!--December 4th, 2020-->

<details>
<summary><font size="+1"><b>December 4th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>We now support reactions from Chinese Hamster Ovary (CHO) for integrated pathway analysis in <a href="https://docs.elucidata.io/Apps/Multi-omic%20Data/IntOmix.html">IntOmix</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Resolved timeout error for opening a folder containing large number of files within a Workspace.</li>
    <li>Resolved issue with Workspace root directory redirection on selection.</li>
  </ul>
</details>

<hr>

<!--November 20th, 2020-->

<details>
<summary><font size="+1"><b>November 20th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Improved OmixWiki UI for better consumption.</li>
    <li>Added the ability to clone Notebooks within Workspaces.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added granular error messages for Notebook functions and CLI jobs.</li>
    <li>Resolved the issue with renaming large data files.</li>
    <li>Resolved the issue with folder breadcrumb in Workspaces.</li>
    <li>Fixed involuntary logout issue.</li>
  </ul>
</details>

<hr>


<!--November 6th, 2020-->

<details>
<summary><font size="+1"><b>November 6th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Data transfer time limit has been extended to 8 hour enabling transfer of 1TB data through <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">CLI</a> at once.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated user interface of Discover and <a href="https://docs.elucidata.io/Apps/Data%20Studio/Data%20Studio.html">Data Studio</a>.</li>
    <li>Added filtering interface to GEO data lake.</li>
    <li>Added search functionality on Discover interface.</li>
    <li>Added highlight and cumulative size feature on multiselection in <a href="https://docs.elucidata.io/Getting%20Started/Workspaces.html">Workspaces</a>.</li>
    <li>Updated collaborators icon to show number of collaborators.</li>
    <li>Resolved inconsistent log<sub>2</sub>FC values for multiple comparisons in <a href="https://docs.elucidata.io/Apps/Multi-omic%20Data/IntOmix.html">IntOmix</a>.</li>
    <li>Resolved sample name descrepancy in concentration plot of <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/QuantFit.html">QuantFit</a>. 
    <li>Fixed table column resizing error on filtering interface.</li>
    <li>Resolved a bug in Polly Docker Domain.</li>
  </ul>
</details>

<hr>

<!--October 23rd, 2020-->

<details>
<summary><font size="+1"><b>October 23rd, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Hosted our first <a href="https://elucidata.io/ugm/">User Group Meeting</a>.</li>
    <li>Introduced our public platform <a href="https://omixwiki.elucidata.io/dashboard">OmixWiki</a>, showcasing top 100 cited COVID-19 publications with end to end omics analysis.</li>
    <li>Released the newest version of <a href="https://github.com/ElucidataInc/ElMaven/releases">El-MAVEN v0.12.0</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated <a href="https://docs.elucidata.io/Getting%20Started/Workspaces.html">Workspaces</a> user interface.</li>
    <li>Added filtering interface to COVID-19 data lake.</li>
    <li>Updated datasets searchability on dataset ID and description.</li>
    <li>Fixed incorrect memory error in <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">CLI</a>.</li>
  </ul>
</details>

<hr>


<!--October 9th, 2020-->

<details>
<summary><font size="+1"><b>October 9th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced the option to make dockers on Polly public by adding public docker domain.</li>
    <li>Welcome screen now displays the username.</li>
    <li>Decreased launch time for applications and notebooks through horizontal pod scaling and buffering.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Fixed landing on Discover after logging in error.</li>
    <li>Fixed priority assignment of automated jobs error.</li>
    <li>Fixed renaming files after upload error.</li>
    <li>Fixed 404 error in Metabolomics Data Lake.</li>
    <li>Integrated documentation to every application.</li>
  </ul>
</details>

<hr>

<!--September 25th, 2020-->

<details>
<summary><font size="+1"><b>September 25th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced Labeled LC-MS Analysis preset for natural abundance correction and visualization for single or dual labeled LC-MS data combined with an interactive, customizable and shareable reporting dashboard.</li>
    <li>Integrated pathway visualization in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html"> Labeled LC-MS Workflow</a>.</li>
    <li>Added dilution factor and protein normalization in the <a href="https://docs.elucidata.io/Apps/Lipidomics%20Data/Lipidomics%20Visualization%20Dashboard.html"> Lipidomics Visualization Dashboard</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added warning message to prevent duplicate folder creation in Workspaces.</li>
    <li>Fixed nested folder creation and notebook renaming error in Workspaces.</li>
    <li>Fixed 503 error in Metabolomics Data Lake.</li>
    <li>Fixed a bug associated with notebooks and shiny apps opening to a blank screen.</li>
    <li>Fixed error occurring in automated jobs.</li>
  </ul>
</details>

<hr>

<!--September 11th, 2020-->

<details>
<summary><font size="+1"><b>September 11th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced <a href="https://docs.elucidata.io/Apps/Data%20Studio/Data%20Studio.html"> Data Studio</a> that brings the tools you need to create, customize, and share your analysis effortlessly with your team across the world.</li>
    <li>Introduced <a href="https://docs.elucidata.io/Apps/Data%20Studio/CCLE%20Correlation%20Analysis.html"> CCLE Correlation Analysis</a> for identification of features correlated with a gene mutation such as mutations in other genes, expression and sample level metadata.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated the version of scanpy to 1.6.0 in single cell docker.</li>
    <li>Fixed a bug in notebook giving error with CLI commands.</li>
  </ul>
</details>

<hr>


<!--August 28th, 2020-->

<details>
<summary><font size="+1"><b>August 28th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced a metabolomics docker equipped with packages for analysis of metabolomics data.</li>
    <li>Added restore functionality to all the <a href="https://docs.elucidata.io/Data%20Lake.html#data-lake-applications"> Data Lake applications</a>.</li>
    <li>Added boxplots for lipids in <a href="https://docs.elucidata.io/Apps/Lipidomics%20Data/Lipidomics%20Visualization%20Dashboard.html"> Lipidomics application</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Updated discoverpy package in all the dockers to the latest version.</li>
    <li>Fixed CellxGene visualization loading for specific datasets.</li>
    <li>Fixed duplicate metabolite generation issue within the <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html"> Dual Mode Data Visualization application</a>.</li>
    <li>Fixed minor UI issues in Workspaces.</li>
    <li>Decreased Workspaces loading time.</li>
  </ul>
</details>

<hr>


<!--August 14th, 2020-->

<details>
<summary><font size="+1"><b>August 14th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced <a href="https://docs.elucidata.io/Getting%20Started/Workspaces.html"> Workspaces</a> on Polly, which is a new and improved version of Polly Projects.</li>
    <li>Added GTEx app to process the filtered datasets from GTEx data lake.</li>
    <li>Added a filtering interface for GTEx data lake that allows filtering of the data on the basis of fields within the curated dataset.</li>
    <li>Integrated <a href="https://docs.elucidata.io/Data%20Lake.html#polly-discover-app"> Discover</a> and <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html"> Dual Mode Visualization</a> for processing and further analysis of transcriptomic and metabolomic and single cell filtered datasets.</li>
    <li>Integrated <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20Notebooks.html"> Notebook</a> to process the filtered datasets.</li>
    <li>Hosted <a href="https://chanzuckerberg.github.io/cellxgene/">CellxGene</a> for processing and visualization of single cell datasets.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Enabled logs access functionality through <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">Polly CLI</a>.</li>
    <li>Added the python package, Discoverpy to all the dockers.</li>
  </ul>
  <p class="Deprecated-button">Deprecated</p>
  <ul>
    <li> The Project Management Dashboard has been deprecated and replaced by Workspaces.</li>
  </ul>
</details>

<hr>



<!--July 31st, 2020-->

<details>
<summary><font size="+1"><b>July 31st, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>Added dot plot for Gene Ontology in the <a href="https://docs.elucidata.io/Data%20Lake.html#polly-discover-app"> Discover</a> application.</li>
    <li>Added an extra layer of security in authentication.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Allowed internal standards and unlabeled data to pass through the <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow</a> to generate output.</li>
    <li>Added Phantasus, Boxplot & Whisker plot along with the bar plot in the <a href="https://docs.elucidata.io/Data%20Lake.html#polly-discover-app"> Discover</a> application.</li>
    <li>Fixed Polly CLI auto login error in notebooks.</li>
    <li>Fixed unresponsive notebook with infinite loading.</li>
  </ul>
</details>

<hr>

<!--July 17th, 2020-->

<details>
<summary><font size="+1"><b>July 17th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>We have released the newest version of <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">Polly CLI v0.1.18</a> enabling you to run a CLI job without the need of "secret" key if the private docker is on Polly.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li><a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow</a> has N and C as indistinguishable isotopes.</li>
    <li>Improved the stability of both Shiny and Desktop Applications.</li>
    <li>Communication within the infrastructure is now through encrypted keys.</li>
    <li>Shiny apps as well as shiny states are encrypted during transit as well as storage.</li>
    <li>Added encryption for the disks running the computations.</li>
    <li>Encrypted buckets containing credentials.</li>
  </ul>
</details>

<hr>

<!--July 3rd, 2020-->
<details>
<summary><font size="+1"><b>July 3rd, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>We have released the newest version of <a href="https://github.com/ElucidataInc/ElMaven/releases/tag/v0.11.0">El-MAVEN v0.11.0.</a></li>
    <li>Polly now provides its own docker repository for easy <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html#docker-management">management of dockers</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added Si as an indistinguishable isotope in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow.</a></li>
    <li>Introduced pre-processing functionalities along with updated selections and heatmap for visualization in <a href="https://docs.elucidata.io/Apps/Lipidomics%20Data/Lipidomics%20Visualization%20Dashboard.html">Lipidomics Visualization Dashboard</a>.</li>
  </ul>
  <p class="Deprecated-button">Deprecated</p>
  <ul>
    <li>Deprecated El-MAVEN FirstView Integration.</li>
  </ul>
</details>

<hr>

<!--June 19th, 2020-->
<details>
<summary><font size="+1"><b>June 19th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
   <li>We now support reactions from <i>Drosophila melanogaster</i> for integrated pathway analysis in IntOmix.</li>
    <li>Introduced <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html#statistical-analysis">pathway enrichment and pathway view</a> feature along with <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html#comparative-analysis">comparative analysis</a> in Dual Mode Data Visualization.</li>
    <li>DEPMAP CCLE (DEPMAP Cancer cell line expression data and dependency scores for genes) repository has been added in <a href="https://docs.elucidata.io/Data%20Lake.html">Data Lake</a>.</li>
    <li>Implemented input file access from the sub-folders of a project for applications.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>The Single Cell Downstream docker is updated with these new packages: rpy2, anndata2ri (Python packages), ExperimentHub (R package).</li>
    <li>Added a GPU instance for <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">Polly CLI</a>.</li>
  </ul> 
</details>

<hr>

<!--June 5th, 2020-->
<details>
<summary><font size="+1"><b>June 5th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
   <li>Introduced visualization of labels in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html#visualization">stacked plot</a> within Labeled LC-MS Workflow.</li>
    <li>Enabled least privilege access for stringent access policies.</li>
    <li>Encryption of data in transit and at rest.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Improved access logs throughout the platform.</li>
    <li>Enhanced security using a secrets management service.</li>
    <li>Implemented regular backups and versioning of data.</li>
  </ul> 
</details>

<hr>

<!--May 22nd, 2020-->
<details>
  <summary><font size="+1"><b>May 22nd, 2020</b></font></summary>
  <br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/CompoundDiscoverer%20QuantFit.html">Polly QuantFit</a> node in <a href="https://www.thermofisher.com/in/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/compound-discoverer-software.html">Compound Discoverer<sup>TM</sup></a> that allows peak picking and absolute quantification on raw data obtained from a Thermo Scientific<sup>TM</sup> Mass Spec instrument.</li>
  </ul>
</details>

<hr>

<!--May 8th, 2020-->
<details>
  <summary><font size="+1"><b>May 8th, 2020</b></font></summary>
  <br>
  <p class="new-button">New</p>
  <ul>
    <li>We now host our desktop application, <a href="https://docs.elucidata.io/Apps/Metabolomic Data/El-MAVEN.html">El-MAVEN on Polly</a>.</li>
    <li><a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MSMS%20Workflow.html#phibeta-tab">Phi calculation</a> feature has been added to Labeled LC-MS/MS Workflow.</li>
  </ul> 
  <p class="update-button">Update</p>
  <ul>
    <li>Changed the optimized color palette in IntOmix from a red-yellow-green scale to a more intuitive red-green scale. All upregulated metabolites or genes are represented by a shade of red and downregulated metabolites or genes as a shade of green.</li>
    <li>Changed the non-optimized color palette in IntOmix from a pink-purple scale to a red-green scale to remove ambiguity.</li>
  </ul> 
</details>  

<hr>

<!--April 24th, 2020-->
<details>
  <summary><font size="+1"><b>April 24th, 2020</b></font></summary>
  <br>
  <p class="new-button">New</p>
  <ul>
    <li>COVID-19 (Transcriptional datasets for SARS viruses, viral infections, and therapeutics for novel coronavirus) repository has been added in <a href="https://docs.elucidata.io/Data%20Lake.html">Data Lake</a>.</li>
    <br>
    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
    <iframe src="https://www.youtube.com/embed/AYgAb5Lbj4g" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
  </ul> 
</details>  

<hr>

<br />

<!--button style-->
<style>
  .update-button {
    background-color: #4C61AF;
    border: 1px solid #364574;
    border-radius: 70px;
    color: #FFFFFF;
    padding: 0px 5px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 12px;
    margin: 4px 2px;
    cursor: default;
  }
  .new-button {
    background-color: #4CAF50;
    border: 1px solid #367437;
    border-radius: 70px;
    color: #FFFFFF;
    padding: 0px 5px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 12px;
    margin: 4px 2px;
    cursor: default;
  }
  .Deprecated-button {
    background-color: #b30000;
    border: 1px solid #b30000;
    border-radius: 70px;
    color: #FFFFFF;
    padding: 0px 5px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 12px;
    margin: 4px 2px;
    cursor: default;
  }
</style>
