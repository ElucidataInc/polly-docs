Polly CLI jobs can now be initiated, managed and have a status-checked for from Polly Python. This lets users run jobs on the Polly cloud infrastructure by scaling computation resources as per need. Users can start and stop tasks and even monitor jobs.

### submit_job()

With this, the job will be submitted to run and Job ID will be created. This Job ID will be needed to check the status and the logs of the submitted job.

```
from polly.auth import Polly;
from polly.jobs import jobs; 
job = jobs()

Polly.auth(AUTH_TOKEN)
job = jobs()

job_file = "<json_job_file>"
workspace_id = <worspace_id>

job.submit_job(workspace_id,job_file)
```

Argument description:-

- workspace_id (str/int): the id of the workspace where the job has to submitted.
- job_file (str) : a json file path which contains the description of a job

Example job file
```
{
  "cpu": "100m",
  "memory": "64Mi",
  "image": "docker/whalesay",
  "tag": "latest",
  "name": "exampleName",
  "command": [
      "cowsay",
      "hello world"
  ]
}
```

### job_cancel()

This function is used to cancel an ongoing job.

```
from polly.auth import Polly;
from polly.jobs import jobs; 
Polly.auth(AUTH_TOKEN)
job = jobs()
job.job_cancel(workspace_id, job_id)
```

Argument description:-

- workspace_id (str/int): the id of the workspace where the job has to submitted.
- job_id (str) : job id to be cancelled


### job_status()

This function is to be used for checking status of a job.

```
from polly.auth import Polly;
from polly.jobs import jobs; 
job = jobs()

Polly.auth(token_s)

job = jobs()
job.job_status(workspace_id, job_id)
```

Argument description:-

- workspace_id (str/int): the id of the workspace where the job has to submitted.
- job_id (str) : job id
 

Checking status of all jobs in the workspace

```
from polly.auth import Polly;
from polly.jobs import jobs; job = jobs()
Polly.auth(AUTH_TOKEN)

job = jobs()
job.job_status(workspace_id)
```

Argument description:-

- workspace_id (str/int): the id of the workspace.

