While creating a docker to be run on Polly, the following must be taken care of.

*   Dockers must be present in either Docker Hub or Amazon ECR. 

*   Soon you will be able to have Dockers directly on Polly.

*   Only self contained dockers can be run on Polly. A self contained docker is one which has the code to get input files as well as upload output files back contained in the docker.

*   Public as well as private dockers are supported. In order to run private dockers, “secret” should be passed as a key in the JSON file. If your private dockers are on Polly itself, you don't require to generate this secret.

*   To get the secret key for the private docker, the following steps need to be followed.

    *   For MacOS, you need to remove the key value pair "credsStore": "osxkeychain" from the config.json file present in the directory `/Users/< username >/.docker`.

    *   You need to be logged in to DockerHub or ECR through the terminal. If not, you will need to log in.

    *   Run the command `sudo polly` on the terminal. 

    *   Select the option miscellaneous followed by create secret for docker. 

    *   Provide the path to the docker config file (the usual path for docker config is `/Users/< username >/.docker/config.json` in Mac and `/home/< username >/.docker/config.json` in Linux). Relative paths are not supported. 

    *   Select the account in which the docker to be run is present. 

    *   Copy the long text string (secret key) output to the JSON file in the key “secret”.

<pre><code>{
  "cpu": 1,
  "memory": "1Gi",
  "image": "docker/whalesay",
  "tag": "latest",
  "Secret": "ewoAImF1dVnphR0ZzWjNWd2RHRTZSVkJKUXlOcFlXMGsiCgkw==",
  "name": "exampleName",
  "command": [
      "cowsay","hello world"
  ]</code></pre>

*   **Passing Environment variables:** Two types of Environment variables can be passed in the json file.

    *   **Normal environment variables** are saved in a database for future references. These can be passed in the parameter **"env”** in the json file.

    *   **Private environment variables** are not saved in any database. These can be used for passing credentials in the json file. These can be passed in the parameter **“secret_env”** in the json file.

**Note:**

*   The value of Environment variables should always be string. For example, the correct way to assign Environment variable is `{“parallel_threads” : “2”}` and **NOT** `{“parallel_threads” : 2}`.

<pre><code>{
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
}</code></pre>

*   [Here](https://gist.github.com/GeorgeSabu/8a3251e263d93b08413ce2c56d8af45d "https://gist.github.com/GeorgeSabu/8a3251e263d93b08413ce2c56d8af45d") is an example gist showing how input data for a job can be taken from and output stored back to Polly Workspaces.

##Polly CLI help

If help is needed for any command, just type `--help` at the end of the command and execute.

![Polly CLI Help](../img/PollyCLI/8.png "Polly CLI Help") <center>**Figure 13.** Polly CLI Help</center>


