# Release Notes
<!--August, 2022-->

<details open>
<summary><font size="+1"><b>August, 2022</b></font></summary>
<br>
 
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
  
Updates:- 

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
