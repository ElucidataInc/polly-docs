#Introduction

##Overview

CRISPR-Cas9 (Clustered Regularly Interspaced Short Palindromic Repeats) is one of the most popular gene-editing techniques to modify gene loci at the desired position in a model organism. The gene-editing CRISPR/cas9 method cuts at specific locations of a  DNA with the help of a small stretch of guide RNAs that takes the Cas9 endonuclease to a specific site. Increasingly the scientific groups are applying this technology to create genetic screens in order to identify mutations that drive treatment resistance to cancer or rapidly access the drug targets.  The aim of the genome-wide CRISPR screening experiment is to screen mutant cells to identify genes involved in a particular phenotype. The beauty of the pooled screen is to integrate many genetic perturbations in one experiment. The CRISPR/cas9 system begins by generating a library of perturbed cells using the library of gRNAs. Lentivirus or Retrovirus is used to deliver gRNAs i.e. integrated to the genome of the cell and serve as a molecular tag. After that cells are separated according to the phenotype of interest and the genes causing the phenotype change can be read out by first isolating the genomic DNA from the population of the cell using PCR followed by sequencing(NGS) across gRNA encoding regions. Then computationally sequencing read count is mapped to a precompiled list of the designed gRNA library. 

**How does CRISPR screening help scientists find gene targets?**

CRISPR-Cas9 (Clustered Regularly Interspaced Short Palindromic Repeats) is one of the most popular gene-editing techniques to modify gene loci at the desired position in a model organism. CRISPR technique is extremely helpful in the identification of genes involved in a particular biological phenomenon. The two main components of CRISPR screening are Cas9 (CRISPR-associated protein 9) and specific guide RNAs that either disrupt host gene or insert gene fragments of interest. Bacteria use this system as a part of their adaptive immune response. Cas9 is an enzyme involved in this mechanism and sgRNA is the complementary sequence of foreign DNA. The popularity of CRISPR-Cas9 is also influenced by its simplicity and cost-effectiveness.

CRISPR Screening provides an interactive way to explore and visualize the quality and analysis results of CRISPR screens. Presently, the tool analyzes the data using​ MAGeCK-VISPR algorithm​ and display the results in an interactive and intuitive manner.

##Scope of the app

*   Swiftly upload large data sets and analyze multiple comparisons in one go.
*   Track the processing status of each comparison.
*   Visualize and interpret the quality check and analysis results effectively.
*   Compare the results of different analyses on the same dashboard.
*   Save and share analyses with other collaborators.
*   Revisit previously performed analyses at any time.

![CRISPR](../img/CRISPR/QuantFit.png) <center>**Figure 1.** QuantFit</center>

#Getting Started

##User Input

CRISPR Screening requires the following three files as input:

**File 1 Description**

**File 2 description**

**File 3 description**

##Steps involved for data processing

*   **Plan a run:** To obtain insights from your data generated from a CRISPR screening experiment, open the CRISPR Screening app on Polly.
*   **Create/Choose a project:** Before running your data in CRISPR Screening, you need to either create or choose an existing project on Polly.
*   **Uploading data:** Polly CRISPR Screening requires two types of input files; fastq.gz files generated from the experiment and the sgRNA library.
*   **Select pipeline and parameters:** You need to select parameters such as *Pipeline*, *Species*, *Assembly*, *Trim*, and *sgRNA length*.
*   **Rename samples:** You can rename the samples if required.
*   **Add comparisons:** After filling all parameters you need to create comparisons and define baseline samples as well as various conditions.
*   **Start processing:** Once the comparisons are set you can start processing.
*   **Log files:** You can also look into the progress of run by analyzing log files.
*   **Go to the dashboard:** Once the run has completed you can click on go to the dashboard and analyze data.
*   **Set base comparison:** Select the base comparison to see the table of genes with the *p*-value and beta score.
*   **Quality control:** You can check the quality of analysis at the sequence level as well as the sample level. In sequence level quality check, you can find graphs of sequencing reads, GC content, Base quality, Mapped reads, Gini index, and zero counts. Similarly at the sample level there are PCA clustering, Read count CDF and normalization graphs.
*   **MAGeCK VISPR:** MAGeCK VISPR output consists of sgRNA Read count plot and any gene of interest can be selected in various samples. Similarly, Hierarchical clustering over normalized beta scores (using Ward variance minimization and Euclidean distance) on 500 genes is also provided. You can also analyze the various comparisons and distribution plots.
*   **Cross-analysis comparison:** In Cross-analysis comparison feature you can select any relevant comparison from previous analysis and compare with the results present analysis.

##Caveats

#Tutorial

##Uplaod Files

##Step 2

##Step 3

##Step 4

##Step 5

#Details about the app

##Detail 1

##Detail 2

##Detail 3

#References

#Frequently Asked Questions (FAQs)