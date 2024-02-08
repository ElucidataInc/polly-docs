Polly CLI can run dockerized jobs on managed Polly infrastructure. Polly infrastructure will scale computational resources with increased usage. All you need to do is submit a job and rest is taken care of by Polly. 

## Create job description JSON file

JSON file is needed to describe the job to be run on Polly. This file should contain the information about the computational resources (machine), docker image, the name of the job and specific commands (if required) to be run after the docker has been run, as keys. Text can be copy pasted from the example below to create the JSON file.

<pre><code>{
 "machineType" : "gp",
 "image": "docker/whalesay",
 "tag": "latest",
 "name": "Single Cell RNA",
 "command": [
     "cowsay","hello world"
 ]
}</code></pre>

**MachineType**

Name of the machine required to run the job needs to be mentioned as per the following table.

| machineType | No. of vCPUs | Memory (RAM) | No. of GPUs |
|-------------|--------------|--------------|-------------|
| ci2xlarge | 16 | 32 GB | - |
| ci3xlarge | 36 | 72 GB | - |
| mi2xlarge | 4 | 32 GB | - |
| mi3xlarge | 8 | 64 GB | - |
| mi4xlarge | 16 | 122 GB | - |
| mi5xlarge | 32 | 250 GB | - |
| mi6xlarge | 64 | 500 GB | - |
| mi7xlarge | 64 | 970 GB | - |
| gpusmall | 4 | 16 GB | 1 |
| gpumedium| 32 | 240 GB | 4 |


If you need a specific machine to be added to the list, please contact us at [polly.support@elucidata.io](mailto:polly.support@elucidata.io).

Depending on the computational power needed, use the key **“machineType”** in the JSON file.

*   **image**: The path to the docker image present in DockerHub or ECR needs to be mentioned in this key.

*   **tag:** Tag of the docker image needs to be mentioned in this key.

*   **name:** Name you want to provide to the job has to be mentioned in this key.

*   **command:** Any commands to be executed after the docker has been run can be mentioned in this key.


## Polly CLI help

If help is needed for any command, just type `--help` at the end of the command and execute.

![Polly CLI Help](../img/PollyCLI/8.png "Polly CLI Help") <center>**Figure 13.** Polly CLI Help</center>


## Some useful gists

*   [Accessing Polly files in and out of a job](https://gist.github.com/GeorgeSabu/8a3251e263d93b08413ce2c56d8af45d)

*   [Running a cluster of jobs with different parameters](https://gist.github.com/GeorgeSabu/e89891da1d86fbaa3afa0655a4ede899)

*   [Bash script to identify when a job finishes](https://gist.github.com/GeorgeSabu/4fbc359fa9ee2bf4d3cb05df3b60db81)
