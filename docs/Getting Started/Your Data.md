#Adding files and notebooks to your workspace

###Using GUI

To add files to your workspaces, click on the *cloud upload icon* present above the middle panel. You can choose the *Upload Files* or *Upload a Notebook* option from within the menu in order to upload the data files and notebooks respectively.

**Note:**

*    If your file size is greater than **100 MB**, use the Polly CLI option to upload to the workspace.

![File Addition](../img/Workspace/7.png) <center>**Figure 8.** File Addition</center>

###Using CLI


You can use PollyCLI as well to import your data directly within Polly. This is a more convenient option for uploading large data (more than **100 MB** in size). The steps for the same are detailed [here]( https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html).


##Deleting items within a workspace

Deleting Multiple Files: Select the files that you want to delete through the checkbox. The selection menu will appear on the top right. Click on *Delete icon*.

![Deleting Workspace](../img/Workspace/21.png) <center>**Figure 22.** Deleting Multiple Files</center>

Confirm by clicking inside the checkbox and then click on *Delete* to delete your selected files.

![Deleting Workspace](../img/Workspace/22.png) <center>**Figure 23.** Deletion Permission</center>

Deleting a single file: If you need to delete a single file, an alternate option can also be used. Click on the kebab menu present at the right end and choose the *Delete* option.

![Deleting Workspace](../img/Workspace/23.png) <center>**Figure 24.** Deleting a single file</center>


##Renaming a file within workspace

Click across the file you want to rename on the kebab menu. Select the *Rename* option.

![Renaming Workspace](../img/Workspace/24.png) <center>**Figure 25.** Renaming Workspace</center>


Provide the new name to the file and click on *Rename* to confirm your changes.

![Renaming Workspace](../img/Workspace/25.png) <center>**Figure 26.** Confirm Renaming</center>


##Accessing and Downloading a file within a workspace

*    Using GUI


To access and download files from your workspaces, Click on the file to look at the file details on the right panel. You can look at the versions of the file along with general details like the file size and modifications date.

**Note:**  

*    If your file size is greater than **100 MB**, use the Polly CLI option to access and download it.

You can download the file by choosing the *Download* option provided at the bottom of the panel.

![Accessing and Downloading a file](../img/Workspace/26.png) <center>**Figure 27.** Accessing and Downloading a file</center>

*    Using CLI


You can use PollyCLI as well to import your data directly within Polly. This is a more convenient option for accessing large data (more than **100 MB** in size). The steps for the same are detailed [here]( https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html).


##Accessing and Restoring an analysis within a workspace

Click on the analysis to look at the file details on the right panel. You can look at your *saved states* as versions of your analysis along with input details. 

![Accessing and Restoring an analysis](../img/Workspace/27.png) <center>**Figure 28.** Accessing and Restoring the saved states</center>

Choose the version through the checkbox and click on the *Restore icon* to restore the saved state.

![Accessing and Restoring an analysis](../img/Workspace/29.png) <center>**Figure 29.** Accessing and Restoring the analysis</center>

You can also choose to look at the state history by clicking on the *History icon*. It gives a detailed description of the parameters and the input files used to save the selected state. Restore after reassuration. 

![Accessing and Restoring an analysis](../img/Workspace/29.5.png)

![Accessing and Restoring an analysis](../img/Workspace/30.png) <center>**Figure 30.** Accessing and Restoring the analysis</center>

For some applications, you can select the algorithm to look at the parameters used for the selected analysis. Click on the *History icon* and then *Retore* the selected algorithm.

![Accessing and Restoring an analysis](../img/Workspace/31.png) <center>**Figure 31.** History icon</center>

![Accessing and Restoring an analysis](../img/Workspace/32.png) <center>**Figure 32.** Restoring algorithm</center>


##Accessing and Restoring a Pollyglot notebook within a workspace

Select the Polly Notebook you want to launch. You can look at the notebook details on the right panel. If you are running the notebook for the first time, the option *Edit and Launch* would appear as a default selection to launch the selected notebook. You are required to select an environment and a machine to run the given notebook, once these selections are done you can launch the notebook.

![Accessing and Restoring a Notebook](../img/Workspace/33.png) <center>**Figure 33.** Accessing and Restoring a Notebook</center>

Select the desired *Docker Environment* from the drop-down menu.

![Accessing and Restoring a Notebook](../img/Workspace/34.png) <center>**Figure 34.** Docker Environment Selection</center>

as well as the *Machine Type* required to run the job. Once the selections are done, click on *Launch*.

![Accessing and Restoring a Notebook](../img/Workspace/35.png) <center>**Figure 35.** Machine Type Selection</center>

For an older notebook: You have two options, you can either launch the notebook directly by the *Launch* button or you can choose to edit it first before launching through the *Edit and Launch* button.

**Note:** 

*    Under the edit option, you can only change the machine type. The docker environment would remain the same as the one selected when you run the notebook for the first time.

![Accessing and Restoring a Notebook](../img/Workspace/36.png) <center>**Figure 36.** Launch for old notebooks</center>


##Starting a New Analysis

Click on the *New Analysis* option provided at the top of the right panel. You will be redirected to the Application interface where you can select an application to launch.

![Starting a New Analysis](../img/Workspace/37.png) <center>**Figure 37.** Starting a New Analysis</center>


##Accessing Documentation

Click on the “?” icon present at the bottom left corner to access documentation.

![Accessing Documentation](../img/Workspace/38.png) <center>**Figure 38.** Accessing Documentation</center>
