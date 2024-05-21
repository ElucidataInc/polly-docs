
Welcome to PWL quick start guide!

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
demo_protein_processing/        # pwl pipeline named "demo_protein_processing"
    │
    ├──  __init__.py
    │
    ├── build
    │   ├── Dockerfile          # For building docker image (must be present)
    │   └── requirements.txt    # dependencies for pipeline (must be present)
    │
    ├── config                  
    │   ├── dev.json            # config for devpolly
    │   ├── test.json           # config for testpolly
    │   └── prod.json           # config for polly
    │
    ├── src                     # Source code
    │   ├── __init__.py
    │   └── main.py 
    │
    └── parameter_schema.json   # Defines pipeline's parameters (must be present)
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
    mkdir pipelines/pwl/<name_your_pipeline>
    cp -r pipelines/pwl/demo_protein_processing/ pipelines/pwl/<name_your_pipeline>/
    ```

3. Go to `build/Dockerfile` and change the pipeline path in the highlighted `COPY` command

    ``` hl_lines="3"
    FROM mithoopolly/workflows-base:python3.9

    COPY pipelines/pwl/demo_protein_processing/build/requirements.txt .

    RUN pip install -r requirements.txt

    ```

4. Change the entrypoint function name in `main.py` to match pipeline name. This is important!

    ``` python hl_lines="2 12"
    @workflow(result_serialization=Serialization.JSON)
    def demo_protein_processing(exp_id: str = "exp1", pre_process: bool = False):
        secret_key = "MY_SECRET_KEY"
        secret_value = Secrets.get(secret_key)
        Logger.info(f"My secret value: {secret_value}")

    ##
    ##
    ##

    if __name__ == "__main__":
        demo_protein_processing("exp1.data", True)
    ```


5. After all the above changes are done, let's push your pipeline

    ``` bash
    git add .
    git commit -m 'First pipeline'
    git push origin <name_of_your_branch>
    ```

6. Go to [circleCI](https://app.circleci.com/pipelines/github/ElucidataInc/pipelines) and approve the hold to deploy your pipeline 


Congrats! You have deployed your first pipeline. Go to [Polly](https://polly.elucidata.io/manage/pipelines) (after circleCI jobs are completed). Click on your pipeline, pass in the parameters and initiate your first run. 
 

Now that you have deployed your first PWL pipeline, let's do in-depth dive on creating your pipelines from scratch. [Check this page](UnderstandingTheSyntax.md).

<br>
<br>
<br>
<br>
<br>
<br>