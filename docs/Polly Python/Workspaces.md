
Polly python enables the users to connect OmixAtlas with Workspaces. Currently, there are two functions to create a new workspaces and listing the existing workspaces. The following library needs to be imported for users to work with workspaces.

```
from polly.workspaces import Workspaces
```

### create_workspace()
Create a new workspace with desired name
```
workspaces.create_workspace("name_of_workspace")
```
### fetch_my_workspaces()
Fetch existing workspaces
```
workspaces.fetch_my_workspaces()
```  
##### upload_to_workspaces()
Upload files or folder to a workspace
```
workspaces.upload_to_workspaces(workspace_id = int, workspace_path = str, local_path = str)
```
### download_to_workspaces()
Download files or folder from a workspace
```
workspaces.download_to_workspaces(workspace_id = int, workspace_path = str)
```
### create_copy()
This function enables user to create a copy of a file or folder contained in a workspace. The copy of the file/folder gets created in the specified destination workspace.

```
from polly.workspaces import Workspaces
workspaces = Workspaces(token)

workspaces.create_copy(source_id, source_path, destination_id, destination_path)
```
`source_id`(int) : workspace id of the source workspace where the file/folder exists.

`source_path` (str): file/folder path on the source workspace to be copied.

`destination_id`(int) : workspace id of the destination workspace where the file/folder is to be copied.

`destination_path` (str):optional parameter to specify the destination path.
