## Single cell RNA-seq workflow on Polly Notebook:

The Polly notebook is a Jupyter Polyglot notebook that supports multiple kernels across different cells of the same notebook. Polly CLI is pre-installed in the docker and can be used here. For the sake of this use case, we will refer to the Single cell dataset and the pipeline recommended [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html "https://satijalab.org/seurat/articles/pbmc3k_tutorial.html").

1\. To begin with the analysis, download the raw data given at the link above and upload it in your workspace through Polly GUI as shown below:

![Gists_Figures](../img/Gists_Figures/SCUpload.png) <center>**Figure .7** Polly workspace</center>

Once the dataset is uploaded, it will be visible in the workspace as shown below:

![Gists_Figures](../img/Gists_Figures/SCfile.png) <center>**Figure .7.2** Uploaded data in workspace</center>

2\. Launch a new Polly notebook in your workspace with the desired docker environment and machine type

![Gists_Figures](../img/Gists_Figures/Singlecellnb.png) <center>**Figure .8** Polly notebook</center>

3\. Read the dataset from your workspace to a notebook instance to make it accessible for analysis as below using Polly CLI in Bash kernel:

![Gists_Figures](../img/Gists_Figures/Gists_fig9.png) <center>**Figure .9** Polly files copy</center>

4\. Uncompress the file read in the notebook instance:

![Gists_Figures](../img/Gists_Figures/Gists_fig10.png) <center>**Figure .10** Untar data</center>

Please note that the uncompressed directory filtered\_gene\_bc\_matrices and its contents will be lost on closing the notebook session. If the uncompressed folder is required to be saved (although not necessary here), it can be done as below:

![Gists_Figures](../img/Gists_Figures/Gists_fig11.png) <center>**Figure .11** Polly files sync</center>

You can also rename the folder while syncing with the workspace. Further, on refreshing the workspace page, you will be able to see the compressed data added along with your notebook and the downloaded data.

![Gists_Figures](../img/Gists_Figures/Gists_fig12.png) <center>**Figure .12** Polly workspace with data</center>

5\. Follow the reference link and perform the scRNA-seq analysis as per the steps described in R. One step is shown below as an example.

![Gists_Figures](../img/Gists_Figures/Gists_fig13.png) <center>**Figure .13** Seurat analysis</center>

6\. In case a particular R library of interest is not available in the pre-built environment, you can easily install it in the notebook instance using the notebook itself or the Terminal from Polly Offerings menu as shown below:

![Gists_Figures](../img/Gists_Figures/Gists_fig14.png) <center>**Figure .14** Polly Offerings</center>

Enter the R console (note the usage of sudo) and install the desired library. Here, as an example, we will install library(pryr). First, we can check that it is missing by calling the library in R as below:

![Gists_Figures](../img/Gists_Figures/Gists_fig15.png) <center>**Figure .15** Missing R library</center>

You can install the library from Bioconductor as below:

![Gists_Figures](../img/Gists_Figures/Gists_fig16.png) <center>**Figure .16** Installing R library</center>

7\. Please remember that the output files generated during the analysis will also need to be saved/copied in the workspace specifically, else they will be lost on exiting the session. The individual files can be copied using the polly files copy command as shown below:

<pre><code>polly files copy -y -s “output.txt“ -d polly://Folder/output.txt</code></pre>

8\. Finally, a recommendation for long running jobs; please run them using the standalone Polly CLI tool to be able to launch them in the backend.

Enjoy analyzing your data with Polly CLI!!
