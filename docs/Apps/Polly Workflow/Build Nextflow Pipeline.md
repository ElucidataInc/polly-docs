## Build Nextflow Pipeline

### Get the workspace ready:

1. Visit https://bitbucket.org/elucidatainc/nextflow/src/k8s-onboarding/ - Can't find link, download the repository folder. Now the files and folders in the repository holds the following details:

    ![](../img/Polly_Workflow/build-nextflow-pipeline-repository.png)

2. Update the data, modules and[  main.nf](http://main.nf/) files as per your requirement - resource to understand this better can be found here -[Nextflow training](https://training.seqera.io/#_simple_rna_seq_pipeline)

    - **Data** - Data folder should be updated with all the relevant fq files
    - **Module** - Its basically a stand-alone module scripts that can be included and shared across multiple workflows. Each module can contain its own `process` or `workflow` definition.


    Each pipeline can be divided into stand alone modules, where the DSL language can be used to define the input, and output structure:
  
    Example: 
    ```    
        process fastqc {
            cpus 2

            publishDir "${params.fastqc_outdir}", mode: 'copy', overwrite: false

            input:
            tuple val(sample_id), file(reads)

            output:
            file("*") 

            script:
            """
            fastqc -t 2 -o "./" -f fastq -q ${reads}

            """  
        }
        
    ```
        
    Understand how modules work via the .nf files under modules folder in the repository

    - [Main.nf](http://main.nf/) - Is the main nextflow script which imports all the modules to run. This allows you to store these components in a separate file(s) so that they can be reused in multiple workflows.

    Sections under the[  main.nf](http://main.nf/) file:

    #### 2.1. Define your params:

    ```
    params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
    params.fastqc_outdir= "FASTQC" 
    params.multiqc_outdir= "MULTIQC"
    ```

    #### 2.2. Define your polly workspace: 

    ```
    params.workspace_id = 9810
    ```

    Import Modules :
    ```
    include {fastqc} from './module/fastqc.nf' params(params)
    include {multiqc} from './module/multiqc.nf' params(params)
    include {pollySync_fastqc} from './module/pollySync_fastqc.nf' params(params)
    include {pollySync_multiqc} from './module/pollySync_multiqc.nf' params(params)
    ```

    #### 2.3. Message: 

    ```
    log.info """

    Nextflow on K8s with DSL2 modules | TUTORIAL
    =============================================
    pipeline reads : ${params.reads}
    fastqc output directory : ${params.fastqc_outdir}
    multiqc output directory : ${params.multiqc_outdir}
    polly workspace id : ${params.workspace_id}
    =============================================
    """
    .stripIndent()
    ```

    #### 2.4. Main script:
    ```
    workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
    fastqc(read_pairs_ch)
    pollySync_fastqc(fastqc.output.collect())
    multiqc(fastqc.output.collect())
    pollySync_multiqc(multiqc.output.collect())
    }
    ```

    #### 2.5. Completion handler
    ```
    workflow.onComplete { 
	println ( workflow.success ? "\nWorkflow completed. Enjoy! \n" : "Oops .. something went wrong" )   
    }
    ```

    #### 3. Create a folder in your workspace, give a relevant name, here we have used demo as a name, upload the downloaded and updated repository files and folder from step  1 and 2 from local to Polly workspace.

    <p align="center">Or</p>

    ##### Using Polly CLI function the local files can be synced to the workspace

    To upload the downloaded and updated repository from step  3 and 4 from local to Polly workspace, by using  the following command:

    ```
    polly files sync --workspace-id  --source  --destination 
    ```
    Workspace ID of the workspace where the data is being synced has to be mentioned in the `--workspace-id` option. Source and destination can be Polly workspace path as well as local path. Workspace path should start with `polly://` followed by the directory path in the workspace where the data is to be synced. Here `polly://` is the root directory for the mentioned workspace.

    ```
    polly files sync --workspace-id 9810 --source ./ --destination polly://Demo
    ```
    







