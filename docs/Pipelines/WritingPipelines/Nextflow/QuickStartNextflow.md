
Welcome to Nextflow quick start guide!

## Setting up the environment

Start by cloning the repository. Assuming you have your ElucidataInc GitHub SSH key setup:
``` bash
git clone git@github.com:ElucidataInc/pipelines.git
```

Create a virtual environment in Python [refer to this doc](https://www.freecodecamp.org/news/how-to-setup-virtual-environments-in-python/) and activate it

``` bash
cd pipelines
# Activate your virtual env here
```

Install some basic requirements

``` bash
pip install -r requirements.txt
```

Install pre-commit hooks for basic formatting checks on code commit

``` bash
pre-commit install
```

<hr>


## Understanding the structure of pipelines repo

Let’s go over the structure of the repository in brief. The following schematic shows some important root level files and folders and their purposes:

``` hl_lines="7 8 9 10 11 12 13 14"
pipelines                       # the repository
    │
    ├── .circleci/              # config for CI/CD
    ├── deployment/             # deployment scripts and utilities
    ├── orchestration/          # utilities for enabling pipeline development
    │
    ├── pipelines/          
    │   │
    │   ├── nextflow/           # All Nextflow pipelines
    │   │   ├── pipeline_1/     
    │   │   └── pipeline_2/     
    │   │
    │   └── pwl/                # All PWL pipelines
    │       └── pipeline_3/     
    ├── ...
    ├── requirements.txt        # dependencies
    ├── ...
    └── scripts/                # common scripts
```

!!! info
    As a pipeline developer, you should only care about the pipelines directory (highlighted above). It will contain both Nextflow and PWL pipelines

<hr>

A Pipeline will follow a specific directory structure. To better grasp this concept, let's explore the directory structure of a demo pipeline.

```
toy/                            # nextflow pipeline named "toy"
    │
    ├──  __init__.py
    │
    ├── build
    │   ├── Dockerfile          # For building docker image (must)
    │   └── environment.yml     # dependencies for pipeline (must)
    │
    ├── config                  
    │   ├── dev.json            # config for devpolly
    │   ├── test.json           # config for testpolly
    │   └── prod.json           # config for polly
    │
    ├── src                     # Source code
    │   ├── main.nf
    │   ├── Makefile
    │   └── nextflow.config   
    │
    └── parameter_schema.json   # Defines pipeline's parameters (must)

```


## Let's create your first pipeline

1. We will start by forking a branch from `#! master`

    ``` bash
    git checkout master
    git checkout -b <add_your_branch_name>_dev
    # Make sure your branch name ends with _dev.
    ```

    The pipelines repository employs a branching strategy. For more details please refer to [this page](../../BranchingStrategy.md).


2. Secondly, instead of creating a pipeline from scratch, let's copy an example pipeline and try playing with it

    ``` bash
    mkdir pipelines/nextflow/<name_your_pipeline>
    cp -r pipelines/nextflow/toy/ pipelines/nextflow/<name_your_pipeline>/
    ```

3. Go to `build/Dockerfile` and change the pipeline path in the highlighted `COPY` command

    ``` hl_lines="4"
    FROM nfcore/base:2.1

    # Install the conda environment
    COPY pipelines/nextflow/toy/build/environment.yml .
    RUN pip3 --no-cache-dir install --upgrade awscli

    CMD ["bash","echo 'ECS_IMAGE_PULL_BEHAVIOR=once' >> /etc/ecs/ecs.config"]
    ```

4. After all the above changes are done, let's push your pipeline

    ``` bash
    git add .
    git commit -m 'First pipeline'
    git push origin <name_of_your_branch>
    ```

6. Go to [circleCI](https://app.circleci.com/pipelines/github/ElucidataInc/pipelines) and approve the hold to deploy your pipeline 


Congrats! You have deployed your first pipeline. Go to [Polly](https://polly.elucidata.io/manage/pipelines) (after circleCI jobs are completed). Click on your pipeline, pass in the parameters and initiate your first run. 
 

<br>
<br>
<br>
<br>
<br>
<br>