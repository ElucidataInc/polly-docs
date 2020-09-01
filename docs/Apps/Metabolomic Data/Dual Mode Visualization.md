#Introduction

##Overview

The Dual Mode Data Visualization (Metabolomics App) allows you to perform downstream analysis on single mode (either positive or negative mode) as well as dual mode (both positive and negative mode) targeted, semi-targeted (without retention time) and untargeted unlabeled metabolomics data along with insightful visualizations. The app provides a variety of normalization methods, scaling options and data visualization functionalities, thereby allowing an efficient analysis of the data to get actionable insights.

##Scope of the App

*   The application supports data with a simple matrix having samples in the columns and metabolites in the rows.
*   It provides different normalization and scaling methods to perform on the data.
*   Performs quality checks for internal standards, metabolites, and samples.
*   Performs statistical analysis using limma and provides interactive visualizations.
*   Provides heatmap visualization along with different algorithms like hierarchical clustering, k-means, correlation etc.
*   Performs comparative analysis for the different cohort comparisons.

![Dual Mode Data Visualization (Metabolomics App) workflow schematic](../../img/DualMode/metab_app_workflow_schematics.png) <center>**Figure 1.** Dual Mode Data Visualization </center>

#Getting Started

##User Input

To process single mode data, the following files are required:

*   El-MAVEN Output File

![El-MAVEN raw intensity file obtained using a compound database](../../img/DualMode/metab_app_single_mode_elmaven.png) <center>**Figure 2.** El-MAVEN raw intensity file obtained using a compound database</center>

*   El-MAVEN Internal Standards File (Optional)

![El-MAVEN raw intensity file of internal standards only](../../img/DualMode/metab_app_int_stds_elmaven.png) <center>**Figure 3.** El-MAVEN raw intensity file of internal standards</center>

*   Cohort Mapping File

![Cohort file](../../img/DualMode/metab_app_single_mode_metadata.png) <center>**Figure 4.** Cohort file</center>

To process dual mode data, the following files are required:

*   El-MAVEN Output Files from both positive and negative mode

![El-MAVEN raw intensity files for positive and negative modes](../../img/DualMode/metab_app_single_mode_elmaven.png)

![El-MAVEN raw intensity files for positive and negative modes](../../img/DualMode/metab_app_dual_mode_neg_elmaven.png) <center>**Figure 5.** El-MAVEN raw intensity files for positive and negative modes</center>

*   El-MAVEN Internal Standards File from both, positive and negative modes (Optional)

![El-MAVEN raw intensity files for positive and negative modes of internal standards only](../../img/DualMode/metab_app_int_stds_elmaven.png)

![El-MAVEN raw intensity files for positive and negative modes of internal standards only](../../img/DualMode/metab_app_dual_mode_neg_stds_elmaven.png) <center>**Figure 6.** El-MAVEN raw intensity files for positive and negative modes of internal standards only</center>

*   Cohort Mapping File

![Sample-Cohort metadata file](../../img/DualMode/metab_app_dual_mode_metadata.png) <center>**Figure 7.** Sample-Cohort metadata file/center></center>

**NOTE:**

*   An already processed .gct file can also serve as input to the app.
*   The internal standard file is optional.

#Caveats
Pathways analysis only works when the data has KEGG Ids within the “compoundId” column.

#Tutorial

Select *Dual Mode Data Visualization (Metabolomics App)* from the dashboard under the *Metabolomics Data* Tab as shown in Figure 8. Create a *New Workspace* or choose from the existing one from the dop-down and provide the *Name of the Session* to be redirected to Dual Mode Data Visualisation (Metabolomics App)'s upload page.

![Polly Dashboard and Workspace selection](../../img/DualMode/Dashboard.png)

![Polly Dashboard and Workspace selection](../../img/DualMode/Selection.png) <center>**Figure 8.** Polly Dashboard and Workspace selection</center>

##Upload Files

The Upload Files interface allows you to upload the input files required for processing through the app which includes the group summary matrix files from El-MAVEN and metadata file, or a .gct file.

![Upload Files Tab](../../img/DualMode/metab_app_uploads_tab.png) <center>**Figure 9.** Upload Files Tab</center>

*   *Download all sample files:* This option would allow you to download the demo files which include the El-MAVEN output in group summary matrix format and cohort mapping file.

*   *Upload GCT file:* Checking this box prompts the app that the input data is in the .gct format.

![Upload Files Tab: *Upload gct file* option is checked](../../img/DualMode/metab_app_uploads_gct_option.png) <center>**Figure 10.** Upload Files Tab: *Upload gct file* option is checked</center>

*   *Select mode:* This drop-down allows you to select the specific mode of the data in case of single mode, meaning whether it is:
    *   Positive, or
    *   Negative

![Upload Files Tab for single-mode data: default view](../../img/DualMode/metab_app_single_mode_uploads.png) <center>**Figure 11.** Upload Files Tab for single-mode data: default view</center>

To upload the output file of internal standards, click on *Upload internal standards file*.

![Upload Files Tab for single-mode data: *Upload internal standards file* option is checked](../../img/DualMode/metab_app_single_mode_stds_option_uploads.png) <center>**Figure 12.** Upload Files Tab for single-mode data: *Upload internal standards file* option is checked</center>

*   *El-MAVEN positive/negative mode file:* This allows you to upload the positive/negative mode El-MAVEN output file (depending on the mode selected) in the .csv peak table format.

*   Internal Standards positive/negative mode file: This allows you to upload the positive/negative mode internal standards El-MAVEN output file (depending on the mode selected) in the .csv peak table format.

*   *Metadata file:* This allows you to upload the cohort mapping file in the .csv format.

*   *Data is dual mode:* Checking this box prompts the app that the input data is from dual-mode (positive and negative modes).

![Upload Files interface for dual-mode data](../../img/DualMode/metab_app_dual_mode_default_upload.png) <center>**Figure 13.** Upload Files interface for dual-mode data</center>

To upload the output file of internal standards, click on *Upload internal standards file*.

![Upload Files Tab for dual-mode data: *Upload internal standards* *file* option is checked](../../img/DualMode/metab_app_dual_mode_stds_option_uploads.png) <center>**Figure 14.** Upload Files Tab for dual-mode data: *Upload internal standards* *file* option is checked</center>

*   *El-MAVEN negative mode file:* This allows you to upload the negative mode El-MAVEN output file in the .csv peak table format.

*   *El-MAVEN positive mode file:* This allows you to upload the positive mode El-MAVEN output file in the .csv peak table format.

*   *Internal Standards negative mode file (optional):* This allows you to upload the negative mode internal standards El-MAVEN output file in the .csv peak table format.

*   *Internal Standards positive mode file (optional):* This allows you to upload the positive mode internal standards El-MAVEN output file in the .csv peak table format.

*   *Metadata file:* This allows you to upload the cohort mapping file in the .csv format.

Click on *Go* to proceed to the next step.

**Note:**

*   The format of the metadata file for dual-mode should be in a specific format.

To make common sample names across the different modes, click on *Replace samples to common sample names*.

![Make common sample names: *Replace samples to common sample names* option is checked](../../img/DualMode/metab_app_common_samplenames.png) <center>**Figure 15.** Make common sample names: *Replace samples to common sample names* option is checked</center>

##Pre-processing

The Pre-processing interface allows you to perform a multitude of functions on the data such as:

![Preprocessing Tab](../../img/DualMode/metab_app_preprocessing_tab.png) <center>**Figure 16.** Preprocessing Tab</center>

*   *Select Internal Standards:* This allows you to select the internal standard(s) from within the El-MAVEN output file when a separate internal standards file is not provided as input.

**Note:**

*   In case, the internal standard(s) are not in the El-MAVEN output file but in the separate internal standards file, they will not show up in the drop down menu. To select the desired internal standards, select them in *Normalize by individual internal standards* option under *Normalize by Internal standards* in the *Perform Normalization > Normalization*.

![Selecting internal standards from the data](../../img/DualMode/metab_app_int_stds_from_elmaven.png) <center>**Figure 17.** Selecting internal standards from the data</center>

*   *Drop Samples:* This allows you to drop/remove certain samples from further analysis which could be blank samples or any samples that didn’t have a good run during MS processing. Samples can be dropped by clicking on *Drop Samples* as shown in Figure 18 after selecting the sample(s) from the drop down menu.

![*Drop Samples* option](../../img/DualMode/metab_app_drop_samples_option.png) <center>**Figure 18.** *Drop Samples* option</center>

*   *Normalize by Internal standards* performs normalization using the internal standards.

![Normalize by internal standards options](../../img/DualMode/metab_app_norm_by_int_stds.png) <center>**Figure 19.** Normalize by internal standards options</center>

*   *Normalize by sum of internal standards* normalizes by the sum of the standards provided.

*   *Normalize by average of internal standards* normalizes by the average of the standards provided.

*   *Normalize by individual internal standards* normalizes by the internal standards selected previously.

*   *Normalize by metabolites* normalizes by any particular metabolite selected.

![Normalize by metabolites option](../../img/DualMode/metab_app_normalise_by_metab_option.png) <center>**Figure 20.** Normalize by metabolites option</center>

*   *Normalize by sum of metabolites* normalizes by the sum of metabolites. Here, the user can select the metabolites from the dropdown option.

*   *Normalize by metadata column* normalizes by any additional column specified in the metadata file. such as cell number etc.

*   *Normalize by control* normalizes by control samples present in the data.

![Normalize by metadata options](../../img/DualMode/metab_app_normalize_by_metadata.png) <center>**Figure 21.** Normalize by metadata options</center>

*   *log2*

*   y + log2(x) [where data is shifted by max value of data plus one]

![Scaling options](../../img/DualMode/metab_app_scaling_options.png) <center>**Figure 22.** Scaling options</center>

**Note:**

*   If internal standard(s) have already been selected in the *Select Internal Standards* option, they would be present in the drop down.

Clicking on *Run* will perform the normalization and scaling based on the parameters selected.

*   *Table:* This displays the data table and visualizations for both pre- and post- normalization.
    *   **Metadata:** This displays the metadata uploaded. This data can be downloaded in the .csv format as shown in Figure 23.
    *   **Metabolite Mapping data:** This displays the metabolite data uploaded. This data can be downloaded in the .csv format as shown in Figure 24.
    *   **Raw data:** This displays the raw El-MAVEN data uploaded. This data can be downloaded in the .gct format as shown in Figure 25.
    *   **Processed data:** This displays the normalized El-MAVEN data based on the parameters selected. This data can be downloaded in .gct format as shown in Figure 26.
    *   **Pre-Processing Results:** This allows you to have a look at the sample distribution with the help of density plot and box-plot before normalization as shown in Figure 27.
    *   **Post-Processing Results:** This allows you to have a look at the sample distribution with the help of the density plot and box-plot after normalization. This provides you with the ability to check the effect of the normalization parameters on the data as shown in Figure 28.

![Metadata table](../../img/DualMode/metab_app_metadata_table.png) <center>**Figure 23.** Metadata table</center> 

![Metabolite Mapping data table](../../img/DualMode/metab_app_metabolite_mapping_data_table.png) <center>**Figure 24.** Metabolite Mapping data table</center> 

![Raw data table](../../img/DualMode/metab_app_raw_data_table.png) <center>**Figure 25.** Raw data table</center>  

![Processed data table](../../img/DualMode/metab_app_normalized_data_table.png) <center>**Figure 26.** Processed data table</center>  

![Pre-Processing Results](../../img/DualMode/metab_app_pre_norm_results.png) <center>**Figure 27.** Pre-Processing Results</center>

![Post-Processing Results](../../img/DualMode/metab_app_post_norm_results.png) <center>**Figure 28.** Post-Processing Results</center>

##Quality Checks

This tab allows you to perform quality checks for the internal standards, metabolites and across samples with the help of interactive visualizations.

**Internal Standards**

It allows you to have a look at the quality of the internal standards used in the data with the help of the different visualizations for any individual as well as for all internal standards.

*   Internal Standards (Individual): You can visualize the quality checks for any internal standard specifically. This allows you to select the internal standard by name, followed by another drop down to select by uniqueId of the feature. It’s also possible to specify the cohort order for the plots. For dual mode data, you can specify the internal standard of the particular mode from the *Select uniqueIds* drop down.

![Internal Standards (Individual) options](../../img/DualMode/metab_app_internal_standards_individual_page.png) <center>**Figure 29.** Internal Standards (Individual) options</center>

![CV Distribution across cohorts](../../img/DualMode/metab_app_indi_stds_cv_distri_across_cohorts.png) <center>**Figure 30.** CV Distribution across cohorts</center>

![CV Distribution across samples](../../img/DualMode/metab_app_indi_stds_cv_distri_samples.png) <center>**Figure 31.** CV Distribution across samples</center>

*   Metabolites: It allows you to have a look at the quality of the metabolites present in the data with the help of the Coefficient of Variation plots

    *   Metabolites CoV Boxplot visualizes the Coefficient of Variation across different cohorts in the data in the form of the boxplot. It’s also possible to specify the cohort order for the plots as shown in Figure 32.
    *   Metabolites CoV Barplot visualizes the Coefficient of Variation as a quality check for any specific metabolite. To use this, select the metabolite followed by the unique id of the feature using the drop downs shown in Figure 33. It’s also possible to specify the cohort order for the plots as shown in FIgure 33.

![Metabolites CoV Boxplot option and CV Distribution across Cohorts boxplot](../../img/DualMode/metab_app_metabolites_cov_boxplot_options.png) <center>**Figure 32.** Metabolites CoV Boxplot option and CV Distribution across Cohorts boxplot</center> 

![Metabolites CoV Barplot](../../img/DualMode/metab_app_metabolites_cov_barplot.png) <center>**Figure 33.** Metabolites CoV Barplot</center>

![CV Distribution across cohorts for selected metabolite](../../img/DualMode/metab_app_metabolites_cv_barplot_cv_across_cohorts.png) <center>**Figure 34.** CV Distribution across cohorts for selected metabolite</center>

![CV Distribution across samples for the selected metabolite](../../img/DualMode/metab_app_metabolites_cv_barplot_cv_across_samples.png) <center>**Figure 35.** CV Distribution across samples for the selected metabolite</center>

**PCA**

This allows you to understand the clustering pattern between biologically grouped and ungrouped samples.

*   PCA (2D) provides PCA visualization in a two-dimensional manner by selecting the PC values for *x-* and *y-* axes. It’s also possible to specify the cohort order for the plots.

![Two-dimensional PCA plot](../../img/DualMode/metab_app_2d_pca.png) <center>**Figure 36.** Two-dimensional PCA plot</center>

*   PCA (3D) provides PCA visualization in a three-dimensional manner by selecting the PC values for *x-*, *y-* and *z-* axes. It’s also possible to specify the cohort order for the plots.

![Three-dimensional PCA plot](../../img/DualMode/metab_app_3d_pca.png) <center>**Figure 37.** Three-dimensional PCA plot</center>

##Statistical Analysis

This interface allows you to perform differential expression analysis with the aim to identify metabolites whose expression differs between any specified cohort conditions. The 'limma' R package is used to identify the differentially expressed metabolites. This method creates a log<sub>2</sub> fold change ratio between the two experimental conditions and an 'adjusted' *p*-value that rates the significance of the difference.

![Statistical Analysis interface](../../img/DualMode/metab_app_stats_analysis_page.png) <center>**Figure 38.** Statistical Analysis interface</center>

The following parameters are available for selection:

*   Select *Cohort A* and *Cohort B:* Default values are filled automatically for a selected cohort condition, which can be changed as per the cohorts of interest.
*   Select *p-val* or *adj. p-val*: Select either *p*-value or adj. *p*-value for significance.
*   *p-val* or *adj. p-val* cutoff: By default, the value is 0.05 but can be changed if required.
*   *log2FC:* Specify the cut-off for log<sub>2</sub> fold change with the help of the slider.

Once the parameters are specified, click on the *Update* button to plot the volcano plot. Based on the parameters specified, a volcano plot is displayed. The volcano plot helps in visualizing metabolites that are significantly dysregulated between two cohorts.

![Volcano plot](../../img/DualMode/metab_app_volcano_plot.png) <center>**Figure 39.** Volcano plot</center>

*Filtered Metabolites Visualization* provides the visualization of cohort-based distribution of the metabolites that are significant based on the parameters specified.

![Filtered Metabolites Visualization](../../img/DualMode/metab_app_stats_analysis_filtered_metabolites_visualization.png) <center>**Figure 40.** Filtered Metabolites Visualization</center>

*Filtered Normalized Table* contains the normalized data of the metabolites that are significant based on the parameters specified.

![Filtered Normalized Table](../../img/DualMode/metab_app_stats_analysis_filtered_normalised_table.png) <center>**Figure 41.** Filtered Normalized Table</center>

*Filtered Differential Expression Table* contains only the metabolites that have significant *p*-values as specified.

![Filtered Differential Expression Table](../../img/DualMode/metab_app_stats_analysis_filtered_diff_exp_table.png) <center>**Figure 42.** Filtered Differential Expression Table</center>

*Differential Expression Table* contains all the differentially expressed metabolites without any filtering.

![Differential Expression Table](../../img/DualMode/metab_app_stats_analysis_diff_exp_table.png) <center>**Figure 43.** Differential Expression Table</center>

*Pathway Enrichment Analysis* performs the pathway enrichment analysis for the significant metabolites based on the parameters specified for the particular cohort comparison. Click on the *Perform Pathway Analysis* button.
As a result, you get Metabolite Set Enrichment Analysis and Pathway Topology Analysis plots that can be downloaded under the *Plot* panel. You can also obtain the tablular representation of the plots by selecting onto the *Table* panel.

![Metabolite Set Enrichment Analysis Plot](../../img/DualMode/metab_app_stats_analysis_msea_plot.png) <center>**Figure 44.** Metabolite Set Enrichment Analysis Plot</center>

![Pathway Topology Analysis Plot](../../img/DualMode/metab_app_stats_analysis_pathway_topology_plot.png) <center>**Figure 45.** Pathway Topology Analysis Plot</center>

![Metabolite Set Enrichment Analysis Table](../../img/DualMode/metab_app_stats_analysis_msea_table.png) <center>**Figure 46.** Metabolite Set Enrichment Analysis Table</center>

![Pathway Topology Analysis Table](../../img/DualMode/metab_app_stats_analysis_topology_table.png) <center>**Figure 47.** Pathway Topology Analysis Table</center>

*Pathway View* plots the pathway view of the metabolites that show up in the Metabolite Set Enrichment Analysis. It maps and renders the metabolite hits on relevant pathway graphs. This enables you to visualize the significant metabolites on pathway graphs of the respective metabolites they belong to. You can select your metabolite of interest from the drop-down and click on *Plot*. This will plot the pathway view of the metabolism selected. You can also download the plot as a .png file by clicking onto the *Download Pathview Plot* button.

![Pathway View Plot](../../img/DualMode/metab_app_stats_analysis_pathview_plot.png) <center>**Figure 48.** Pathway View Plot</center>

##Visualization

This interface allows you to visualize the cohort-based distribution of a specific metabolite or a group of metabolites on the basis on its normalized intensity values.

![Visualization tab options](../../img/DualMode/metab_app_viz_tab_options.png) <center>**Figure 49.** Visualization tab options</center>

*   *Enter metabolite:* Select the metabolite(s) of interest from the drop down option.
*   *Select uniqueIds:* You can specifically select the metabolic feature of interest for the metabolite from the drop down option.
*   *Select order of cohort:* You can also specify the particular order of the cohort to visualize the bar plot.

Once the parameters are selected, click on *Load Plots* to plot the bar plot for the metabolite.

![Cohort-wise bar plot with the normalized intensity of selected metabolite](../../img/DualMode/metab_app_barplot_viz.png) <center>**Figure 50.** Cohort-wise bar plot with the normalized intensity of selected metabolite</center>

##IntOmix Input

This tab allows you to generate the input for [IntOmix](https://docs.elucidata.io/Apps/Multi-omic%20Data/IntOmix.html) where you can visualize the significantly altered metabolic network modules between any two experimental conditions.

![IntOmix Input Tab options](../../img/DualMode/metab_app_intomix_input_tab.png) <center>**Figure 51.** IntOmix Input Tab options</center>

Specify two or more cohorts from the *Select cohorts* drop down option for which you want to generate the IntOmix input. Once the required cohorts are selected, click on *Generate* to generate the IntOmix input.

![IntOmix Table for the cohort conditions specified](../../img/DualMode/metab_app_intomix_table.png) <center>**Figure 52.** IntOmix Table for the cohort conditions specified</center>

**NOTE:**

*   At least two cohorts are required to create the input file.

##Heatmap

This tab allows you to produce a heatmap of the processed data, so that you can observe the level of expression in a visual form. Click on *Load Heatmap* button to generate the heatmap.

![Heatmap of the processed data](../../img/DualMode/metab_app_heatmap.png) <center>**Figure 53.** Heatmap</center>

##Comparative Analysis

This tab allows you to perform comparative analysis between a set of cohorts in your data. As a result of which you can visualize the UpSet plot of the unique and overlapping metabolites for the selected cohort comparisons. Further, you can also perform pathway analysis on the metabolites for the set intersections of interest.

* *Comparison Parameters* tab allows you to select the cohorts of interest for which you would want to get the set intersections. You can select the cohorts from the *Select cohorts* drop-down and click on *Run* button. Further, you can also specify the *p*-value cut-off and log<sub>2</sub>FC threshold.

![Comparison Parameters Tab options](../../img/DualMode/metab_app_comparative_analysis_uploads_tab.png) <center>**Figure 54.** Comparison Parameters Tab options</center>

You will get a table as a result of the parameters specified which will have the significant metabolites for the different cohort comparisons along with their corresponding *p*-values and log<sub>2</sub>FC values. You can also download this table as a .CSV file.

![Comparison Table](../../img/DualMode/metab_app_comparative_analysis_comparison_table.png) <center>**Figure 55.** Comparison Table</center>

* *UpSet Plot* tab allows you to visualize the set intersections for the cohort comparisons selected where every comparison consists of the significant metabolites associated with the same. You can select the cohort comparisons of interest from the *Select Cohort Comparison* drop-down which represents all the possible comparisons for the cohorts specified in the previous tab. Click on *Plot* to get the UpSet plot for the specified comparisons.

![UpSet Plot options](../../img/DualMode/metab_app_comparative_analysis_upset_plot_options.png) <center>**Figure 56.** UpSet Plot options</center>

![UpSet Plot](../../img/DualMode/metab_app_comparative_analysis_upset_plot.png) <center>**Figure 57.** UpSet Plot</center>

Along with the plot, you can also get all the constituent metabolites for the respective comparisons in a tabular format that can be downloaded as a .CSV file.

![UpSet Plot table](../../img/DualMode/metab_app_comparative_analysis_upset_table.png) <center>**Figure 58.** UpSet Plot Table</center>

* *Pathway Enrichment Analysis* tab allows you to perform the pathway enrichment analysis for the significant metabolites that show up based on the parameters specified in the *Comparison Parameters* tab for the particular set of cohort comparison. Click on the *Perform Pathway Analysis* button. As a result, you get Metabolite Set Enrichment Analysis and Pathway Topology Analysis plots that can be downloaded under the *Plot* panel. You can also obtain the tablular representation of the plots by selecting onto the *Table* panel.

![Metabolite Set Enrichment Analysis Plot](../../img/DualMode/metab_app_comparative_analysis_msea_plot.png) <center>**Figure 59.** Metabolite Set Enrichment Analysis Plot</center>

![Pathway Topology Analysis Plot](../../img/DualMode/metab_app_comparative_analysis_pathway_topology_plot.png) <center>**Figure 60.** Pathway Topology Analysis Plot</center>

![Metabolite Set Enrichment Analysis Table](../../img/DualMode/metab_app_comparative_analysis_msea_table.png) <center>**Figure 61.** Metabolite Set Enrichment Analysis Table</center>

![Pathway Topology Analysis Table](../../img/DualMode/metab_app_comparative_analysis_topology_table.png) <center>**Figure 62.** Pathway Topology Analysis Table</center>

*Pathway View* plots the pathway view of the metabolites that show up in the Metabolite Set Enrichment Analysis. It maps and renders the metabolite hits on relevant pathway graphs. This enables you to visualize the significant metabolites on pathway graphs of the respective metabolisms they belong to. You can select your metabolite of interest from the drop-down and click on *Plot*. This will plot the pathway view of the metabolite selected. You can also download the plot as a .png file by clicking onto the *Download Pathview Plot* button.

![Pathway View options](../../img/DualMode/metab_app_comparative_analysis_pathway_view_options.png) <center>**Figure 63.** Pathway View options</center>

![Pathway View Plot](../../img/DualMode/metab_app_comparative_analysis_pathway_view_plot.png) <center>**Figure 64.** Pathway View Plot</center>

##Data for MetaboAnalyst

This tab allows you to generate the input for MetaboAnalyst that enables statistical, functional and integrative analysis of metabolomics data by providing a variety of modules for different functionalities.

![Data for MetaboAnalyst Tab options](../../img/DualMode/metab_app_metaboanalyst_page.png) <center>**Figure 65.** *Data for MetaboAnalyst* Tab options</center>

Specify the cohort comparison of interest by sleecting the cohorts from the *Select cohort A* and *Select cohort B* drop downs. Click on *Download data for metaboanalyst* to download the .csv file. Clicking on *Go to metaboanalyst* will redirect you to MetaboAnalyst’s homepage.

**NOTE:**

*   Normalized intensity data (normalization method chosen in pre-processing tab) is downloaded.
