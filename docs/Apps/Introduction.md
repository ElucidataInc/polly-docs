#Overview

The applications on Polly are built to process and visualize experimental data ranging from mass-spec based omics to sequencing-based omics, from dual-mode data visualization to the analysis of CRISPR screening. Despite the variety and specificity of the apps present on Polly, there are a couple of features common across all apps. These features are built to help you get started with minimal effort, analyze and share results and processed data with your collaborators with ease.

##Demo Data

To help you get started even if you do not have your own data, every app on Polly has demo data uploaded as shown in Figure 1, Figure 2 and Figure 3. Moreover the demo data also serves as a reference point for the input files and their format required for each app.

![Demo Data for FirstView](../img/AppIntroduction/01_Application_feature_demodata.png) <center>**Figure 1.** Demo Data for Dual Mode Data Visualization</center>


![Demo data for MetScape](../img/AppIntroduction/02_Application_feature_demodata.png) <center>**Figure 2.** Demo data for Proteomics Workflow</center>


![Demo data for Labeled LC-MS app](../img/AppIntroduction/DemoDataLLCMS.png)<center>**Figure 3.** Demo data for Labeled LC-MS app</center>

##Upload Data

All Polly <!-- ([El-MAVEN Phi Relative LC-MS](../El-MAVEN Phi Relative LC-MS), [FirstView](../FirstView), [MetScape](../MetScape), [QuantFit](../QuantFit), [CRISPR Screening](../CRISPR Screening) and [IntOmix](../IntOmix)) --> <!-- ([Dual Mode Data Visualization](../Dual Mode Data Visualization), [Labeled LC-MS Workflow](../Labeled LC-MS  Workflow), [Labeled LC-MS/MS Workflow](../Labeled LC-MS/MS Workflow), [Untargeted Pipeline](../Untargeted Pipeline), [High Throughput Drug Screening](../High Throughput Drug Screening), [RNA Seq Workflow](../RNA Seq Workflow), [Proteomics Workflow](../Proteomics Workflow) and [Lipidomics Visualization Dashboard](../Lipidomics Visualization Dashboard)) --> applications provide you the option to upload data from local storage as shown in Figure 4 and Figure 5.

![Upload data from local storage for MetScape](../img/AppIntroduction/04_Application_feature_uploaddata.png)<center>**Figure 4.** Upload data from local storage for Dual Mode Data Visualization</center>

![Upload data from local storage for Labeled LC-MS app](../img/AppIntroduction/UploadDataLLCMSLocal.png)<center>**Figure 5.** Upload data from local storage for Labeled LC-MS app</center>

Some of the Polly applications also provide the ability to upload input files from Polly workspace by using the option *Import from Polly* as shown in Figure 6. This simplifies data processing as well as makes Polly the platform where biological data can be stored and processed conveniently. 

![Upload data from workspace for Labeled LC-MS app](../img/AppIntroduction/UploadDataLLCMSPolly.png)<center>**Figure 6.** Upload data from workspace for Labeled LC-MS app</center>

##Restore Analysis

All Polly <!-- ([El-MAVEN Phi Relative LC-MS](../El-MAVEN Phi Relative LC-MS), [FirstView](../FirstView), [MetScape](../MetScape), [QuantFit](../QuantFit), [CRISPR Screening](../CRISPR Screening) and [IntOmix](../IntOmix)) --> applications <!-- ([Dual Mode Data Visualization](../Dual Mode Data Visualization), [Labeled LC-MS Workflow](../Labeled LC-MS  Workflow), [Labeled LC-MS/MS Workflow](../Labeled LC-MS/MS Workflow), [Untargeted Pipeline](../Untargeted Pipeline), [High Throughput Drug Screening](../High Throughput Drug Screening), [RNA Seq Workflow](../RNA Seq Workflow), [Proteomics Workflow](../Proteomics Workflow) and [Lipidomics Visualization Dashboard](../Lipidomics Visualization Dashboard)) --> contain the restore functionality that allows any analysis to be restored to the last step. Analyses can be restored by navigating to the desired workspace. Click on the specific analysis and select the algorithm you want to restore. Clicking on *Restore* from the right panel will take you back to the application with the same data. Restore helps laboratories and organizations with standardization across labs and improves reproducibility.


Some of the Polly applications also provide the ability to save at that exact moment of your choice using the *Save the State* option.

![Restore for Shiny applications](../img/AppIntroduction/LLCMSRestore.png) <center>**Figure 7.** Saving the states in Shiny applications</center>

For the applications providing the *Save the  state* option, you have an additional step to choose which state to restore. Multiple states, if saved are depicted as versions which makes it extremely easy to visually identify the state of interest as shown in Figure 8.


##Download Plots & Output

All applications on Polly allow you to download processed data and plots when displayed so that if needed you can verify from a third party about exactly how the data is being processed at a specific step. The tabular data can be downloaded as .csv files whereas the plots can be downloaded as .png, .jpeg or .svg files as shown in Figure 8, FIgure 9 and Figure 10.

![Download Pathway Dashboard as an image or the processed data as a .csv file](../img/AppIntroduction/DownloadMS1.png)<center>**Figure 8.** Download Pathway Dashboard as an image or the processed data as a .csv file</center>


![Download PCA Plot as an image](../img/AppIntroduction/DownloadMS2.png) <center>**Figure 9.** Download PCA Plot as an image</center>


![Download the processed data as a .csv file](../img/AppIntroduction/DownloadLLCMS.png) <center>**Figure 10.** Download the processed data as a .csv file</center>

As with upload, data can be downloaded to local storage for Angular and Shiny apps as well as Polly workspace for Shiny apps.
