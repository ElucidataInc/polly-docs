### Reading Nextflow Output

Note the workflow ID - go to the workspace you had entered before - you will see a folder created with the workflow ID there which will have the html and trace files.

![](../img/Polly_Workflow/reading-nextflow-output-workspaces.png)

Here you’ll find four different files:

1. Dag file - Execution process flow of the pipeline can be viewed, observed if the flow is

<Details>
    <Summary>You can see an example below:
    </Summary>
    ![](../img/Polly_Workflow/reading-nextflow-output-dag-file.png)
</Details>
2. Report : Nextflow creates an HTML execution report: a single document which includes many useful metrics about a workflow execution. The report is organised in the three main sections: Summary, Resources and Tasks.
<ul>
    <li>Summary :The Summary section reports the execution status, the launch command, overall execution time and some other workflow metadata

<Details>
<Summary>You can see an example below</Summary>
![](../img/Polly_Workflow/reading-nextflow-output-report.png)
 </Details>
    <li>Resource Usage: The Resources section plots the distribution of resource usage for each workflow process using the interactive plotly.js plotting library.
        <ul>
            <li>Plots are shown for CPU, memory, job duration and disk I/O. They have two (or three) tabs with the raw values and a percentage representation showing what proportion of the requested resources were used. These plots are very helpful to check that task resources are used efficiently.
    
<Details>
<Summary>You can see an example below</Summary>

![](../img/Polly_Workflow/reading-nextflow-output-resource-usage.png.png)
</Details>
        </ul>
        <li>Tasks: The Tasks section lists all executed tasks, reporting for each of them the status, the actual command script, and many other metrics.

<Details>
<Summary>You can see an example below</Summary>
![](../img/Polly_Workflow/reading-nextflow-output-tasks.png)
</Details>
            </ul>
            
3. Timeline: Nextflow can render an HTML timeline for all processes executed in your pipeline. 

<Details>
<Summary>You can see an example below</Summary>
![](../img/Polly_Workflow/reading-nextflow-output-timeline.png)

        Each bar represents a process run in the pipeline execution. The bar length represents the task duration time (wall-time). The colored area in each bar represents the real execution time. The grey area to the left of the colored area represents the task scheduling wait time. The grey area to the right of the colored area represents the task termination time (clean-up and file un-staging). The numbers on the x-axis represent the time in absolute units eg. minutes, hours, etc.
</Details>


