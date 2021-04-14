Polly CLI can run dockerized jobs on managed Polly infrastructure. Polly infrastructure will scale computational resources with increased usage. All you need to do is submit a job and rest is taken care of by Polly. 

##Create job description JSON file

JSON file is needed to describe the job to be run on Polly. This file should contain the information about the computational resources (machine), docker image, the name of the job and specific commands (if required) to be run after the docker has been run, as keys. Text can be copy pasted from the example below to create the JSON file.

<pre><code>{
 "machineType" : "gp",
 "cpu": 1,
 "memory": "1Gi",
 "image": "docker/whalesay",
 "tag": "latest",
 "name": "Single Cell RNA",
 "command": [
     "cowsay","hello world"
 ]
}</code></pre>

**machineType**

Name of the machine required to run the job needs to be mentioned as per the following table.

| machineType | No. of vCPUs | Memory (RAM) | No. of GPUs |
|-------------|--------------|--------------|-------------|
| gp | 4 | 16 GB | - |
| ci2xlarge | 16 | 32 GB | - |
| ci3xlarge | 36 | 72 GB | - |
| mi2xlarge | 4 | 32 GB | - |
| mi3xlarge | 8 | 64 GB | - |
| mi4xlarge | 16 | 122 GB | - |
| gpusmall | 16 | 61 GB | 1 |

More machines (including some with GPUs) will be added soon. If you need a specific machine to be added to the list, please contact us at [polly@elucidata.io](mailto:polly@elucidata.io).

If computational power required is less than 2 vCPUs and 8 GB RAM, use the keys **“cpu”** and **“memory”** in the JSON file instead of the key **“machineType”**. If all 3 keys are present, **“machineType”** takes priority and the machine will be assigned accordingly. In the example JSON (image) mentioned above, machine selected will be **“gp”** with 4 vCPUs and 16 GB RAM and **NOT** 1vCPU and 1 GB RAM.

*   **cpu:** Mention the number of CPUs needed here. For smaller jobs, just a part of the CPU can also be chosen. For example, if 0.1 vCPUs are required for the job, the number of CPUs can be mentioned as **“100m”**. If more than 2 CPUs are required for the job, use the key **“machineType”** to choose the relevant machine instead of **“cpu”** and **“memory”**.

*   **memory:** RAM required needs to be mentioned in text (eg - “1Gi” or “500 Mi”) in this key. If memory needed is more than 8 GB, use the key **“machineType”** to choose the relevant machine instead of **“cpu”** and **“memory”**.

*   **image**: The path to the docker image present in DockerHub or ECR needs to be mentioned in this key.

*   **tag:** Tag of the docker image needs to be mentioned in this key.

*   **name:** Name you want to provide to the job has to be mentioned in this key.

*   **command:** Any commands to be executed after the docker has been run can be mentioned in this key.
