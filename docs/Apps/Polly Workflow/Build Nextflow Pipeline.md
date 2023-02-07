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
    
### Execute-in-terminal/Polly-CLI

1.  Install/update Polly CLI [Polly CLI - Polly Documentation](https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html "https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html") or open up a notebook on polly with the smallest machine possible (as it already has Polly CLI pre-installed)

    <details>
    <summary>Quick installation guide</summary>
    <br>
    <h3>Installation</h3>
    <h4><b>Dependencies Required for Polly CLI</b></h4>
    The following dependencies are required to be installed before installing Polly CLI:

    -   [Node and npm](https://www.npmjs.com/get-npm "https://www.npmjs.com/get-npm"):

        -   Linux: For installation on Linux, follow the steps mentioned [here](https://www.digitalocean.com/community/tutorials/how-to-install-node-js-on-ubuntu-18-04 "https://www.digitalocean.com/community/tutorials/how-to-install-node-js-on-ubuntu-18-04").

        -   Mac: For installation on Mac, follow the steps mentioned [here](https://www.digitalocean.com/community/tutorials/how-to-install-node-js-and-create-a-local-development-environment-on-macos "https://www.digitalocean.com/community/tutorials/how-to-install-node-js-and-create-a-local-development-environment-on-macos").

    <h4><b>Commands to install</b></h4>
    To install Polly CLI, run the following commands on Terminal / Command prompt:

    - Linux: 

        ```
        sudo npm install -g @elucidatainc/pollycli
        ```

    - Mac: 
        ```
        npm install -g @elucidatainc/pollycli
        ```

    </details>

2. Open your terminal or polly notebook code using bash: 

    <details>
    <summary>Currently a small step has to be done for configuring your polly CLI. Please copy paste these codes on a terminal.</summary> 

    ```
    echo 'const axios = require("axios");
    const chalk = require("chalk");
    const pollyEnv = require("./env.json");
    const pollymsg = require("./message");
    const { getHeaders } = require("./pollyheaders");

    const getWorkflowClient = async () => {
        const{ headers } = await getHeaders();
        const workflowClient = axios.create({
            baseURL: `${pollyEnv.computeApi}/workflow`,
            headers
        });

        return workflowClient;
    }

    export const submitWorkflow = async (workspace_id, config) => {
        try {
            const workflowClient = await getWorkflowClient();
            let { pipeline, main_script, container } = config;

            if(!pipeline.endsWith("/")) {
                pipeline = pipeline + "/";
            }

            const body = {
                workspace_id,
                pipeline,
            }

            if(!!container) {
                body.container = container;
            }

            if(!!main_script) {
                body.main_script = main_script
            }

            const res = await workflowClient.post("/", body);
            pollymsg.pollySuccess(`${res.data.message}, run_id = ${res.data.id}`);
            console.log(chalk.bold(`status: polly workflows status --run-id=${res.data.id}`));
            console.log(chalk.bold(`logs  : polly workflows logs --run-id=${res.data.id}`));
            console.log(chalk.bold(`cancel: polly workflows cancel --run-id=${res.data.id}`));
        } catch(e) {
            if(e.response && e.response.data) {
                pollymsg.pollyError(e.response.data.detail);
            }
            pollymsg.pollyError(e.message);
        }
    }

    export const deleteWorkflow = async (run_id) => {
        try {
            const workflowClient = await getWorkflowClient();
            const res = await workflowClient.delete(`/${run_id}`);
            pollymsg.pollySuccess("workflow execution terminated");
        } catch(e) {
            if(e.response && e.response.data) {
                pollymsg.pollyError(e.response.data.detail);
            }
            pollymsg.pollyError(e.message);
        }
    }

    export const getWorkflow = async (run_id) => {
        try {
            const workflowClient = await getWorkflowClient();
            const res = await workflowClient.get(`/${run_id}`);
            if(res.data.status === "failed") {
                console.log(chalk.bold.red(`workflow ${res.data.status}`));
            } else {
                console.log(chalk.bold.green(`workflow ${res.data.status}`))
            }
        } catch(e) {
            if(e.response && e.response.data) {
                pollymsg.pollyError(e.response.data.detail);
            }
            pollymsg.pollyError(e.message);
        }
    }

    export const getWorkflowLogs = async (run_id) => {
        try {
            const workflowClient = await getWorkflowClient();
            const res = await workflowClient.get(`/${run_id}/logs`);
            if(res.data.logs.length === 0) {
                console.log("logs not yet generated");
                return;
            }

            console.log(chalk.bold.italic("logs:"));
            console.log(res.data.logs);
        } catch(e) {
            if(e.response && e.response.data) {
                pollymsg.pollyError(e.response.data.detail);
            }
            pollymsg.pollyError(e.message);
        }
    }'>./workflows.js
    sudo mv ./workflows.js /usr/lib/node_modules/@elucidatainc/pollycli/src/workflows.js
    cat /usr/lib/node_modules/@elucidatainc/pollycli/src/workflows.js
    ```
    </details>

3. Now open up your terminal and type:

    ```
    polly workflows submit \
    --workspace-id <******> \
    --pipeline <name of the folder> \
    --main-script main.nf \
    --container <docker image>
    ```

    Eg:
    ```
    polly workflows submit \
    --workspace-id 9810 \
    --pipeline Demo \
    --main-script main.nf \
    --container docker.polly.elucidata.io/elucidata/ic_exp:nf
    ```
- Run the command
4. You will get 3 commands to see the status or check the logs or to cancel the workflow.

- Just copy paste and run each of these command from the output and you’ll be able to check the status, log or cancel the run.

![](../img/Polly_Workflow/execute-in-terminal-or-polly-cli-logs.png)

![](../img/Polly_Workflow/execute-in-terminal-or-polly-cli-logs-2.png)

#### [Reading the output has been explained in detail at the end](https://docs.elucidata.io/Apps/Polly%20Workflow/Nextflow%20Output.html "#Reading-Nextflow-Output")







