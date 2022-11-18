This will enable users to generate reports, link reports to a dataset in OmixAtlas and fetch reports linked with dataset in an OmixAtlas

### link_report()
Org admins can now link a report (html or pdf file format) present in a workspace with the specified datasets in an OmixAtlas. Once a report is linked to a dataset in OmixAtlas, it can be fetched both from front-end and polly-python.

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas(token)

omixatlas.link_report(repo_key, dataset_id, workspace_id, workspace_path, access_key)
```
`repo_key`(str): repo_name or repo_id of the repository to be linked

`dataset_id`(str): to which the report should be linked

`workspace_id`(int): where the report is located

`workspace_path`(str): specific folder path and file name of the report which should be linked

`access_key`(str): “private" or "public" depending access type to be granted. If public, then anyone with a Polly account with the link will be able to see the report. If private, then only the individuals who have access to the workspace where reports is stored will be able to see them.

This function returns a success message along with the link which can be used to view the report.

### fetch_linked_reports()

This function will enable users to fetch the list of reports linked to a given dataset in an OmixAtlas.

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas(token)

omixatlas.fetch_linked_reports(repo_key, dataset_id)
```
`repo_key`(str): repo_name or repo_id of the repository 

`dataset_id`(str): for which the reports should be fetched

This function returns a dataframe containing information on who added the report, when it was added and the link.

### generate_report()

This is a MVP release to minimise time taken by users to determine relevance of a dataset in their research, we’re enabling auto-generation of reports. These reports will contain dataset and sample level metadata along with some graphs showing how the samples are distributed. It will help them figure out which cohorts could be of interest in a given dataset. This template can be modified and we’re open to user’s feedback.

In future, we’ll enable users to make custom template so that they can support needs of an Enterprise OmixAtlas as well.

This report is available on the user’s local path as well as an option to upload the report to workspaces is given.

```
from polly.omixatlas import OmixAtlas
omixatlas = OmixAtlas(token)

omixatlas.generate_report(repo_key, dataset_id, workspace_id, workspace_path)
```

`repo_key`(str): repo_name/repo_id for which the report is to be generated

`dataset_id`(str): dataset_id for which the report is to be generated.

`workspace_id`(int): workspace_id to where the report is to be uploaded.

`workspace_path`(str) (Optional Parameter): workspace_path to upload the report to.

### delete_linked_report()
This function can be used by Org Admin to delete the file in workspaces with the specified dataset in OmixAtlas.

```
from polly import Omixatlas 
omixatlas = Omixatlas()
omixatlas.delete_linked_report(repo_key: str, dataset_id: str, report_id: str)
```
Argument description:-

repo_key (str): repo_name/repo_id of the repository which is linked.

dataset_id (str): dataset_id of the dataset to be linked.

report_id (str): report id associated with the report in workspaces that is to be deleted. This id can be found when invoking the fetch_linked_report() function.



