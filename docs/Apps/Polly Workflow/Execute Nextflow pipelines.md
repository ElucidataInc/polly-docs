### Execute Nextflow Pipelines

Pre-requisite:
1.  Have access to the repository folder - [https://bitbucket.org/elucidatainc/nextflow/src/k8s-onboarding/ - Can't find link](https://bitbucket.org/elucidatainc/nextflow/src/k8s-onboarding/) if not write request for access from Indraneel.

2.  Pipeline structured in nf pipeline, basically the [main.nf](http://main.nf/ "http://main.nf") and modules should be ready to use

3.  All the relevant Data files

#### **Section 1 : Get the workspace ready:**

1.  Visit [https://bitbucket.org/elucidatainc/nextflow/src/k8s-onboarding/ - Can't find link](https://bitbucket.org/elucidatainc/nextflow/src/k8s-onboarding/), download the repository folder. Now the files and folders in the repository holds the following details:

    ![](../img/Polly_Workflow/build-nextflow-pipeline-repository.png)
2. Update the data, modules and [main.nf](http://main.nf/ "http://main.nf") files as per your requirement - resource to understand this better can be found here - [Nextflow training](https://training.seqera.io/#_simple_rna_seq_pipeline)

    <Details>
    <Summary>
    Understand the contents of each folder in repository
    </Summary>

    -  **Data** - Data folder should be updated with all the relevant fq files

    -   **Module** - Its basically a stand-alone module scripts that can be included and shared across multiple workflows. Each module can contain its own `process` or `workflow` definition.

    -   [**Main.nf**](http://main.nf/ "http://Main.nf")- Is the main nextflow script which imports all the modules to run. This allows you to store these components in a separate file(s) so that they can be re-used in multiple workflows.
    </Details>


    Here under Params section,  you need to update the workspace name in which the repository is updated and uploaded.

    **Define your polly workspace:**

    ```
    params.workspace_id = 9810
    ```

    <Details>
    <Summary>Sections under the main.tf file:
    </Summary>
    2.1. Define your params:
    
    ```
    params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
    params.fastqc_outdir= "FASTQC" 
    params.multiqc_outdir= "MULTIQC"
    ```
    2.2. Define your polly workspace:

    ```
    params.workspace_id = 9810
    ```

    2.3. Message

    2.4. Main Script:

    ```
    workflow {
        read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
        fastqc(read_pairs_ch)
        pollySync_fastqc(fastqc.output.collect())
        multiqc(fastqc.output.collect())
        pollySync_multiqc(multiqc.output.collect())
    }
    ```

    2.5. Completion handler
    
    ```
    workflow.onComplete { 
        println ( workflow.success ? "\nWorkflow completed. Enjoy! \n" : "Oops .. something went wrong" )
    }

    ```


    </Details>

3. Create a folder in your workspace, give a relevant name, here we have used demo as a name, upload the downloaded and updated repository folder from step  3 and 4 from local to Polly workspace,

    <p align="center">Or</p>

    <Details>
    <Summary>Using Polly CLI function the local files can be synced to the workspace
    </Summary>
    To upload the downloaded and updated repository from step  3 and 4 from local to Polly workspace, by using  the following command:

    ```
    polly files sync --workspace-id  --source  --destination 
    ```

    Workspace ID of the workspace where the data is being synced has to be mentioned in the `--workspace-id` option. Source and destination can be Polly workspace path as well as local path. Workspace path should start with `polly://` followed by the directory path in the workspace where the data is to be synced. Here `polly://` is the root directory for the mentioned workspace.

    ```
    polly files sync --workspace-id 9810 --source ./ --destination polly://Demo
    ```
        
    </Details>

#### **Section 2 : [Same as above](https://elucidatainc.atlassian.net/wiki/spaces/~62c251a8ce5a604dbfb37a52/pages/3867050006/Nextflow+-+User+Documentation# "#")**
