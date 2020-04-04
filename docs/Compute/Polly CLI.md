#Polly CLI

##What is Polly CLI?

The Polly CLI (Command Line Interface) is an open source tool that enables you to interact with Polly services using commands in your command-line shell. Polly CLI lets you upload data and run jobs on the Polly cloud infrastructure by scaling computation resources as per need. You can start and stop jobs, monitor them and view logs. 


##Required System Configurations
Polly CLI can work on any Unix based system (Linux and Mac distributions). It does not work on Windows. We will be releasing a Windows version soon. It can be used on local computers as well as cloud instances and servers. 

There are no specific machine configurations required for Polly CLI. It can work on a system with as low as 512 MB RAM and 1 CPU.


###How to install? 

####Dependencies Required for Polly CLI
The following dependencies are required to be installed before installing Polly CLI:

* Node and npm :  

    * Linux : For installation on Linux, follow the steps mentioned [here](https://www.digitalocean.com/community/tutorials/how-to-install-node-js-on-ubuntu-18-04).

    * Mac : For installation on Mac, follow the steps mentioned [here](https://www.digitalocean.com/community/tutorials/how-to-install-node-js-and-create-a-local-development-environment-on-macos).


####Commands to install

To install Polly CLI, run the following commands on Terminal / Command prompt :

* Linux :
```bash
sudo npm install -g @elucidatainc/pollycli
```

* Mac : 
```bash
npm install -g @elucidatainc/pollycli
```

####Commands to uninstall

To uninstall Polly CLI, run the following commands on Terminal / Command prompt :

* Linux : 
```bash
sudo npm install -g @elucidatainc/pollycli
```

* Mac : 
```bash
npm install -g @elucidatainc/pollycli
```

* Note:
    * “sudo” might have to be used before every command while accessing Polly CLI on cloud instance or server.


###Logging in and out of Polly CLI

####Log in
Open the terminal on the system and execute the following command to log in.
```bash
polly login
```

Enter the Polly Username and Password when prompted.

![Polly Login](./img/1.png "Polly Login") <center>**Figure 1:** Polly Login</center>

Once logged in, you will stay logged in the system and won’t need to log in again even if a new terminal is opened or the system is restarted. You will only need to log in again if you manually log out from the system.


####Log out
Execute the following command to log out
```bash
polly logout
```

* Note:
    * Logging out will not stop a running job.


###Create and View Polly Workspaces

####What are Polly Workspaces?
Polly Workspaces are online workspaces that contain data, analyses, code, logs etc for a specific project or experiment. Data is stored and Analysis is performed within a user chosen workspace. More details about workspaces is mentioned here. **< Link to workspaces documentation in the previous section >**


####Creating a new Workspace

To create a new workspace, use the following command.

```bash
polly workspaces create
```

You will be asked to name the Workspace and provide a description. Once the workspace is created, the workspace name and ID will be shown on the screen. This workspace ID will be needed while creating a JSON file for running jobs.

![Create Workspace](./img/2.png "Create Workspace") <center>**Figure 2:** Create Workspace</center>

####Workspaces

To view all the existing Workspaces with access, use the following command.

```bash
polly workspaces list
```

Users will be asked to select which Workspaces to list:
* all: On selecting this, all the Workspaces will be listed.
* latest or oldest: On selecting these the user will be asked to enter the number of workspaces as shown in the image below.

![List Workspaces](./img/3.png "List Workspaces") <center>**Figure 3:** List Workspaces</center>

Workspace ID will be required for transferring data and running jobs.   

###Data Transfer

Polly CLI can be used to transfer large data to and from Polly Workspaces. Upto 5 TBs of data can be transferred in one go. 


####Listing the files and folders in a Polly Workspace Directory
Files and folders within any Polly Workspace can be listed using the following command.

```bash
polly files list --workspace-id <workspaceid> --workspace-path <path/to/workspace/directory>
```

Workspace id can be obtained by using the command “polly workspaces list” as explained in the previous section. The path to the workspace directory has to start with “polly://”. Eg - The following command will list all the files and folders within the folder ABC in the workspace 1234.

```bash
polly files list --workspace-id 1234 --workspace-path polly://ABC/
```

* Note:
    * This command only shows files and folders just one layer within the directory mentioned (just like the “ls” command on terminal).


####Manually sync data to and from Polly

Polly CLI can be used to sync the data between a Polly Workspace and a local directory. Data can be synced manually both ways using the following command.

```bash
polly files sync --workspace-id <workspaceid> --source <path/source/directory> --destination <path/to/destination/directory>
```
Workspace ID of the workspace where the data is being synced has to be mentioned in the --workspace-id option. Source and destination can be Polly workspace path as well as local path. Workspace path should start with “polly://” followed by the directory path in the workspace where the data is to be synced. Here ”polly://” is the root directory for the mentioned workspace.

The following command will sync data from Polly workspace to current local directory.

```bash
polly files sync --workspace-id 1234 --source polly://directory1/ --destination ./
```

The following command will sync data from current local directory to Polly Workspace directory.

```bash
polly files sync --workspace-id 1234 --source ./ --destination polly://directory1/
```

* Note:
    * Only files that have been changed or added new will get transferred using the sync command. The files that remained unchanged after the last sync will not get transferred. This command can only be used for folders or directories (not for individual files). To transfer just a single file to or from Polly, use the “copy” command mentioned in the next section.


###Copy files to and from Polly

Files can be copied to and from a Polly Workspace using the following command.

```bash
polly files copy --workspace-id <workspaceid> --source <path/to/source/file> --destination <path/to/destination/file>
```

This command will copy an individual file from source to destination. The transfer can be from or to Polly Workspace depending on the source and destination defined. Workspace path should start with “polly://” followed by the directory structure within the Workspace.

* Note:
    * Paths need to be enclosed in double quotes (“ ”) if there are spaces or special characters in the path.


###Running Dockerized Jobs

Polly CLI can run dockerized jobs on managed Polly infrastructure. Polly infrastructure will scale computational resources with increased usage. All a user needs to do is submit a job and rest is taken care of by Polly. 


###Create job description JSON file

JSON file is needed to describe the job to be run on Polly. This file should contain the information about the computational resources (machine), docker image, the name of the job and specific commands (if required) to be run after the docker has been run, as keys. Text can be copy pasted from the example below to create the JSON file.

```json
{
 "machineType" : "gp",
 "cpu": 1,
 "memory": "1Gi",
 "image": "docker/whalesay",
 "tag": "latest",
 "name": "Single Cell RNA",
 "command": [
     "cowsay","hello world"
 ]
}
```

machineType
Name of the machine required to run the job needs to be mentioned as per the following table.

| machineType | No. of vCPUs | Memory (RAM) | No. of GPUs |
|-------------|--------------|--------------|-------------|
| gp | 4 | 16 GB | - |
| ci2xlarge | 16 | 32 GB | - |
| ci3xlarge | 36 | 72 GB | - |
| mi2xlarge | 4 | 32 GB | - |
| mi3xlarge | 8 | 64 GB | - |
| mi4xlarge | 16 | 122 GB | - |

More machines (including some with GPUs) will be added soon. If you need a specific machine to be added to the list, please contact us at [polly@elucidata.io](mailto:polly@elucidata.io).

If computational power required is less than 2 vCPUs and 8 GB RAM, use the keys **“cpu”** and **“memory”** in the json file instead of the key **“machineType”**. If all 3 keys are present, **“machineType”** takes priority and the machine will be assigned accordingly. In the example json (image) mentioned above, machine selected will be **“gp”** with 4 vCPUs and 16 GB RAM and NOT 1vCPU and 1 GB RAM.

* **cpu :** Mention the number of CPUs needed here. For smaller jobs, just a part of the CPU can also be chosen. For example, if 0.1 vCPUs are required for the job, the number of CPUs can be mentioned as **“100m”**. If more than 2 CPUs are required for the job, use the key **“machinType”** to choose the relevant machine instead of **“cpu”** and **“memory”**.

* **memory :** RAM required needs to be mentioned in text (eg - “1Gi” or “500 Mi”) in this key. If memory needed is more than 8 GB, use the key **“machinType”** to choose the relevant machine instead of **“cpu”** and **“memory”**.

* **image** : The path to the docker image present in DockerHub or ECR needs to be mentioned in this key.

* **tag :** Tag of the docker image needs to be mentioned in this key.

* **name :** Name that user wants to provide to the job has to be mentioned in this key.

* **command :** Any commands to be executed after the docker has been run can be mentioned in this key.


### Docker Building Guidelines

While creating a docker to be run on Polly, the following must be taken care of.

* Dockers must be present in either Docker Hub or Amazon ECR. 

* Soon users will be able to have Dockers directly on Polly.

* Only self contained dockers can be run on Polly. A self contained docker is one which has the code to get input files as well as upload output files back contained in the docker.

* Public as well as Private dockers are supported. In order to run Private dockers, “secret” should be passed as a key in the json file. To get the secret key for the private docker, the following steps need to be followed.

    * For MacOS, users need to remove the key value pair "credsStore": "osxkeychain" from the config.json file present in the directory “/Users/<username>/.docker”.

    * The user needs to be logged in to DockerHub or ECR through the terminal. If not, the user will need to log in.

    * Run the command “sudo polly” on the Terminal. 

    * Select the option “miscellaneous” followed by “create secret for docker”. 

    * Provide the path to the docker config file (the usual path for docker config is /Users/<username>/.docker/config.json in Mac and /home/<username>/.docker/config.json in Linux).
    Note - Relative paths are not supported. 

    * Select the account in which the docker to be run is present. 

    * Copy the long text string (secret key) output to the json file in the key “secret”.
```json
{
  "cpu": 1,
  "memory": "1Gi",
  "image": "docker/whalesay",
  "tag": "latest",
  "Secret": "ewoAImF1dVnphR0ZzWjNWd2RHRTZSVkJKUXlOcFlXMGsiCgkw==",
  "name": "exampleName",
  "command": [
      "cowsay","hello world"
  ]
}
```

* Passing Environment variables : Two types of Environment variables can be passed in the json file.

    * Normal environment variables are saved in a database for future references. These can be passed in the parameter “env” in the json file.

    * Private environment variables are not saved in any database. These can be used for passing credentials in the json file. These can be passed in the parameter “secret_env” in the json file.

* Note:
    * The value of Environment variables should always be string. For example, the correct way to assign Environment variable is {“parallel_threads” : “2”} and NOT {“parallel_threads” : 2}

```json
{
 "cpu": "100m",
 "memory": "64Mi",
 "image": "your_docker",
 "tag": "latest",
 "env": {
   "ENV1": "ENV_VALUE1",
   "ENV2": "ENV_VALUE2"
 },
 "secret_env": {
   "SECRET_ENV1": "SECRET_ENV_VALUE1",
   "SECRET_ENV2": "SECRET_ENV_VALUE2"
 },
 "name": "docker running"
}
```

### Execute Job
To execute the job, execute the following command
```bash
polly jobs submit
```

On executing this command, you will be asked to enter the id of the workspace where the job should be run and the path to the job description JSON file. With this, the job will be submitted to run and Job ID will be created. This Job ID will be needed to check the status and the logs of the submitted job.

**Note :** You do not need to create a new Workspace for running a job. You can simply list the older Workspaces and run a job in an already created Workspace.

![Submit Jobs](./img/4.png "Submit Jobs") <center>**Figure 4:** Submit Jobs</center>

Monitor Job status

Get job status 
* The following command can be used to view the status of a particular job.
```bash
polly jobs status --workspace-id <workspace id> --job-id <job id>
```

![Single Job Status](./img/5.png "Single Job Status") <center>**Figure 5:** Single Job Status</center>

* The following command can be used to view the statuses of all the jobs in a workspace.
```bash
polly jobs status --workspace-id <workspace id>
```

A prompt to enter job id will appear which when kept blank gets all the job statuses in a workspaces.

![All Job Statuses in a Workspace](./img/6.png "All Job Statuses in a Workspace") <center>**Figure 6:** All Job Statuses in a Workspace</center>

### Get job logs
To view the logs of any job, use the following command.
```bash
polly jobs logs --workspace-id <workspace id> --job-id <job id>
```

This will give the logs for the job. In case the job is still running, it will give the logs generated till that instant.

![Job Logs](./img/7.png "Job Logs") <center>**Figure 7:** Job Logs</center>

### Polly CLI help

If help is needed for any command, just type --help at the end of the command and execute.

![Polly CLI Help](./img/8.png "Polly CLI Help") <center>**Figure 8:** Polly CLI Help</center>
