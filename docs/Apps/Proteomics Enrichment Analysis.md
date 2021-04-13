Proteomics Enrichment Analysis

# Introduction

## Overview

Oncology is in the midst of a paradigm shift towards personalized medicine with the increased wave of next-generation sequence data from labs and clinics. To the modern biologist, the genomic and transcriptomic features have the potential to function as a key guide for selecting patient treatments. Clinical trials remain costly, however, and the need for reliable biomarkers for cancer drug sensitivity is as great as ever.

Proteomics is the scientific discipline, which mainly analyzes the presence and the abundance of proteins in a given sample. Generally, the main components are *expression analysis* as well as *data standardization* and *data publication*. In clinical applications, [protein expression](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/protein-expression) analysis usually aims at the detection of proteins that are differentially expressed between two known groups of patients (e.g., “diseased” vs. “healthy”) and which may be interesting biomarker candidates. To this end, common statistical methods are used. The resulting biomarker candidates must be validated using independent samples and independent technology.

## Scope of the App

-   Annotate the column metadata with a gene expression.

    

-   Filter samples from the expression data and metadata based on selection.

    

-   Analyze the correlation between the expression status of the selected gene against the other features.

    

-   Performs univariate analysis like limma t-test along with multivariate analysis like PCA analysis.

    

-   Measuring transcription factors regulating differentially expressed genes which further associates it to PPIs or Protein-Protein interactions.

    

-   Understand which genes are most suitable for suppressing or exciting a process.

    

-   Performs the enrichment analysis on the gene set relayed.

    

# Tutorial

Select *Proteomics Enrichment Analysis* from the dashboard under the *Studio Presets* tab.

![Polly Dashboard](../img/Notebooks/1.GettingStartedSelectApp.png "Polly Dashboard") <center>**Figure 1.** Polly Dashboard</center>

*Select an existing Workspace* from the drop-down and provide the *Name of the Session* to be redirected to *Transcriptomics Enrichment Analysis* Preset's upload page.

![App Selection](../img/Notebooks/2.SelectWorkspace.png "App Selection") <center>**Figure 1.** App Selection</center>

# Data Curation

## Upload Files

The first component is the *Upload component* which allows you to upload the input files required for processing through the app. It includes the raw expression, sample metadata, and mutation files from Proteomics Repository.

![Select the Dataset you want to explore](../img/Notebooks/3. Global Inputs.png "Select the Dataset you want to explore") <center>**Figure 1.** Select the Dataset you want to explore</center>

You will need to provide the following:

1.  ***Repo ID***: It is the ID of the Repository on Polly that you are referring to that contains your dataset. You will find this ID by either asking the database manager the repo ID or by selecting said dataset in the database and opening it with Jupyter Notebooks. While the notebook is loading, pay attention to the URL link and you will see a substring ‘repo\_id’ followed by a set of numbers. Remember to do this quickly as once the notebook finishes loading, the URL will change to that of the notebook.

    

2.  ***Dataset ID***: That’s the ID of your Dataset. You will find a column by this name when you select a repository.


![DatasetID](../img/Notebooks/4. DatasetID.png "DatasetID") <center>**Figure 1.** DatasetID</center>
    

3.  ***Dataset ID Column***: This field shall be removed in the near future, but for the time being, It is advisable to always set the field to ‘dataset\_id’.

    

Once all the components have been entered, click *Save Global Outputs* (you cannot alter them subsequently if it is successfully accepted).

To start your analysis, click on *Data Picker* and click *Run Task*.

![Execute the Component](../img/Notebooks/5. DataPicker.png "Execute the Component") <center>**Figure 1.** Execute the Component</center>

**Note:** When you launch this session again from your workspace, although you will see the above layout as is and proceed to click *Run Task* in the hope to run the analysis again, the clicking action will not be registered and you will be required to click on *Data Picker* to activate the *Run Task* button action.

The Processing Step will go through 4 phases and will generate three output files:

![Generated Output](../img/Notebooks/8. ThreeOutputFiles.png "Generated Output") <center>**Figure 1.** Generated Output</center>

-   **Column metadata**: Table for the metadata information associated with each sample. Rows represent samples and columns represent various attributes of each sample.

    

-   **Data matrix**: The rows represent different proteins while the columns represent the various samples from whom data was collected. The entries are expressions of each protein.

    

-   **Row metadata**: Table for the metadata information associated with each protein. Rows represent proteins and columns represent various attributes of each Protein.

    

# Data Exploration

## PCA

This component allows you to simplify the complexity of high-dimensional data while retaining the trends and patterns in it. It projects the data onto a lower dimension with an objective to find the best summary of the data using a limited number of principal components that help in understanding the clustering pattern between biologically grouped and ungrouped samples.

![PCA](../img/Notebooks/9. PCA.png "PCA") <center>**Figure 1.** PCA</center>

-   **Cohort Column:** Dropdown to select one of the metadata columns.

    

-   **Top N Variants:** The top N variable entities will be used for PCA calculation. Define the number in this box. The default number used is 1500.

    

It generates two outputs:

-   **PCA Plot:** A plot is created where the samples are labeled based on the cohort selected in the metadata column. When you hover over the points, sample ID and percentage of variance explained by each principal component are displayed along with the cohort.

    

-   **PCA Score:** Table of the first 10 PC values and metadata columns.

![PCA outputs](../img/Notebooks/10. PCAScores.png "PCA outputs") <center>**Figure 1.** PCA outputs</center>

    

## Box Plot

This component generated a *Gene/Metabolite expression Boxplot*. It is a visualization component that uses a five-number summary (minimum, first quartile, median, third quartile, and maximum) to display the distribution of data based on quartiles. Outliers are also displayed as individual points.

![Boxplot](../img/Notebooks/11. BoxPlot.png "Boxplot") <center>**Figure 1.** Boxplot</center>

-   **Select gene:** Dropdown for picking the gene/metabolite/feature for which the boxplot is to be made.

    

-   **Select metadata column:** Dropdown to select one of the metadata columns.

    

After executing the component, an interactive *Boxplot* will be generated as an output for the selected gene/metabolite.

![Boxplot parameters](../img/Notebooks/12. BoxPlotOutput.png "Boxplot parameters") <center>**Figure 1.** Boxplot parameters</center>

## Differential Expression

This component allows the search for differentially expressed (DE) genes, that is, genes that show differences in expression level between conditions or in other ways are associated with given predictors or responses.

![Differential expression parameters](../img/Notebooks/13. DifferentialExpressionParameters.png "Differential expression parameters") <center>**Figure 1.** Differential expression parameters</center>

-   **Cohort Column:** Dropdown to select one of the metadata columns.

    

-   **Cohort A:** Dropdown to select a cohort from the metadata column selected.

    

-   **Cohort B:** Dropdown to select another cohort from the metadata column selected.

    

-   **Normalization:** Perform log2normalization if data is not normalized.

    

-   **Algorithm:** You can select any one of the two algorithms - *limma* or *Unpaired t-test*.

    

*Limma* is an R package for the analysis of gene expression microarray data, especially the use of linear models for analyzing designed experiments and the assessment of differential expression. Limma provides the ability to analyze comparisons between many RNA targets simultaneously in arbitrary complicated designed experiments.

An *unpaired t-test* (also known as an independent t-test) is a statistical procedure that compares the averages/means of two independent or unrelated groups to determine if there is a significant difference between the two.

-   **P-Value Correction:** You can select the *Benjamini-Hochberg* method to correct the *p*\-value for False Discovery Rate or the *Bonferroni* method to correct the *p*\-value for Type I errors.

    

-   **P-Value Metric:** Plot and calculate significance using the selected metric. The *p-value* is the value returned by the algorithm while the *Adjusted p-value* is the corrected value after applying one of the correction methods above.

    

-   **P-value threshold:** You can select the appropriate threshold for the selected *p*\-value metric. *p-values* lower than this threshold will be marked as significant.

    

-   **Absolute Log2FC Threshold:** You can select the appropriate fold change threshold. Log2fold change values higher than this will be marked as significant.

    

Once all the parameters are selected, execute the component by clicking on *Run Task*. It will generate two outputs:

-   **Differential Expression:** Table with Differential Expression results with *p*\-value and fold change.

    

-   **Volcano Plot:** Based on the parameters specified, a volcano plot is displayed. The volcano plot helps in visualizing metabolites that are significantly dysregulated between two cohorts.

![Differential expression outputs](../img/Notebooks/14. VolcanoPlot.png "Differential expression outputs") <center>**Figure 1.** Differential expression outputs</center>


Note: In order to get the visualization, please stay in the *Analyze Data* tab and click on the drop-down on the upper-right corner. Please do NOT open the *Visualization Tab* and expect to see the visualizations there.

## Bar Plot

The *Bar plot* is able to give you a count of the different genes (*Select Gene*) based on a particular way of separating out the samples (*Select Cohort*).

![Bar Plot Input Parameters](../img/Notebooks/15. BarplotInputs.png "Bar Plot Input Parameters") <center>**Figure 1.** Bar Plot Input Parameters</center>

-   **Cohort:** Dropdown to select one of the metadata columns.

    

-   **Gene**: Select one of the Genes to get a count on the basis of each cohort type

    

After executing the component, an interactive *Bar Plot* will be generated as an output.

![Bar Plot output](../img/Notebooks/16. BarplotOutputs.png "Bar Plot output") <center>**Figure 1.** Bar Plot output</center>

# Global Pathway Analysis

## X2K analysis (eXpression2Kinases)

An X2K analysis involves measuring transcription factors regulating differentially expressed genes which further associates it to PPIs or Protein-Protein interactions thereby creating a subnetwork. A Kinase Enrichment analysis is done on the nodes of this subnetwork. 

The X2K analysis is done after the differential expression is carried out. The differentially expressed data is used as an input for X2K analysis. Here, the differential expression is performed where significant genes (*p*\-value < 0.05) are selected. These significant genes are ordered on the basis of their log2FC value. You can select the top 'n' of the ordered values based on the up and downregulation of genes.

You can choose from the following databases to perform the analysis.

-   ChEA 2015

    

-   ENCODE 2015

    

-   ChEA & ENCODE Consensus

    

-   Transfac and Jaspar

    

-   ChEA 2016

    

-   ARCHS4 TFs Coexp

    

-   CREEDS

    

-   Enrichr Submissions TF-Gene co-occurrence

    
![X2K Input Parameters](../img/Notebooks/17. X2KInputParameters.png "X2K Input Parameters") <center>**Figure 1.** X2K Input Parameters</center>

On selecting all the options along with the target database, click on the option *Run Task*. 

![X2K Outputs](../img/Notebooks/18. X2KOutputs.png "X2K Outputs") <center>**Figure 1.** X2K Outputs</center>

## Gene Ontology

Gene Ontology allows us to find out the set of genes that influences a particular process.

One of the previous outputs: Differential Expressions, allows us to understand the effect of a particular stimulus (a drug, in this case) on genes. Some genes may get excited (higher expression or Upregulated) or may get repressed (lower expression or Downregulated). Apart from the altered expression due to the stimulus, we also have expression data for the genes under normal conditions. So, the Differential Expression is able to capture the log-fold change in expression of each gene.

Now, genes play a role in different avenues of an organism. You can select the organism from under the *Select organism name*. The two types that are available are:

-   Homo sapiens,

    

-   Mus musculus.

    

The roles are found under *Ontology Selection*. Under Ontology Selection, we have the following choices:

-   enrichKEGG,

    

-   Molecular Function,

    

-   Biological Process (things like hair-growth), and

    

-   Cellular Components (effects on different cellular components)

![Gene Ontology Input](../img/Notebooks/19. GeneOntologyInput.png "Gene Ontology Input") <center>**Figure 1.** Gene Ontology Input</center>    

*Top Enrichment Pathways* are the set of the top 'n' genes (based on an internal scoring) of the overall effect that gene played in a Process: the top gene was the one whose shutting-down would be most detrimental to a process (maybe because it just so happens that the functions controlled by that gene do not overlap with any other genes functions).

Upon clicking *Run Task*, we get the following output:

## EnrichR (Enrichment analysis)

EnrichR performs the enrichment analysis on the gene set relayed, either in the Differential Expression tab or X2K analysis tab. Enrichment analysis is a computational method for inferring knowledge about an input gene set by comparing it to annotated gene sets representing prior biological knowledge. Enrichment analysis checks whether an input set of genes significantly overlaps with annotated gene sets.

To perform the enrichment analysis, select the desired EnrichR database from [here](https://maayanlab.cloud/Enrichr/#stats) and enter its name in the field.

![EnrichR Inputs](../img/Notebooks/21. EnrichRInputs.png "EnrichR Inputs") <center>**Figure 1.** EnrichR Inputs</center>

## Gene Set Enrichment Analysis

Here, we do the same process as in the case of *Gene Ontology*, except we focus on a specific subset of the gene set (there are, after all, over 20k of them). This list is the same one that was entered in the *EnrichR* process ([source](https://maayanlab.cloud/Enrichr/#stats)).

![Gene Set Enrichment Analysis Inputs](../img/Notebooks/23. GeneSetEnrichmentAnalysisInputs.png "Gene Set Enrichment Analysis Inputs") <center>**Figure 1.** Gene Set Enrichment Analysis Inputs</center>

# Dashboard

You can collect your findings (plots) for different modules using different filters into what is called a [Dashboard](https://docs.elucidata.io/Apps/Data%20Studio/Data%20Studio.html#dashboard). The Visualization Dashboard provides an at-a-glance view of the selected visualization charts. The dashboard is customizable and can be organized in the most effective way to help you understand complex relationships in your data and can be used to create engaging and easy-to-understand reports.

You can add a plot to the Dashboard using the *Add to Dashboard* button.

![Add to Dashboard](../img/Notebooks/25. Dashboard.png "Add to Dashboard") <center>**Figure 1.** Add to Dashboard</center>

Some of the salient features of Dashboard include:

-   Organizing your plots in the order you wish

    

-   Adding plots with different filters

    

-   The generated reports are interactive and can be shared with the collaborators