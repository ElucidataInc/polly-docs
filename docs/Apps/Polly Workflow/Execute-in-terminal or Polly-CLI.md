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

- ![](../img/Polly_Workflow/execute-in-terminal-or-polly-cli-logs-2.png)

#### [Reading the output has been explained in detail at the end](https://elucidatainc.atlassian.net/wiki/spaces/~62c251a8ce5a604dbfb37a52/pages/3867050006/Nextflow+-+User+Documentation#Reading-Nextflow-Output "#Reading-Nextflow-Output")

---