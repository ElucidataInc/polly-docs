class Workspaces:
    """
    This class contains functions to interact with workspaces on Polly. Users can create a workspace, fetch list\
 of workspaces, upload data to workspace and download data from workspace. To get started, users need to \
initialize a object that can use all function and methods of Workspaces class.
    Args:
        token (str): token copy from polly

    Usage:
        # to initialize a obj
        workspaces = Workspaces(token)
    
    """
    def create_workspace(self, name: str, description=None):
        """
        This function create workspace on Polly.
        Args:
              name (str): name of the workspace
              description (str): general information about workspace
        Returns:
              Dictionary (dict) : It will return an object like this
                        {
                        'id': 9999,
                        'name': 'rrrrr',
                        'active': True,
                        'description': 'for docu',
                        'created_time': '2022-03-16 11:08:47.127260',
                        'last_modified': '2022-03-16 11:08:47.127260',
                        'creator': 1127,
                        'project_property': {
                            'type': 'workspaces',
                            'labels': ''
                        },
                        'organisation': 1
                        }

        """

    def fetch_my_workspaces(self):
        """
        This function fetch workspaces from Polly.
        Args:
              None
        Returns:
              it will return a table with attributes
        
        """


    def create_copy(
        self, source_id: int, source_path: str, destination_id: int, destination_path=""
    ) -> None:
        """
        Function to create a copy of files/folders existing in a workspace into another workspace.
        Args:
              source_id (int): workspace id of the source workspace where the file/folder exists
              source_path (str) : file/folder path on the source workspace to be copied
              destination_id (int) : workspace id of the destination workspace where the file/folder is to be copied
              destination_path (str) : optional parameter to specify the destination path
       
        Raises:
              InvalidParameterException: when the parameter like source id is invalid
              InvalidPathException: when the source path is invalid
        
        """


    def upload_to_workspaces(
        self, workspace_id: int, workspace_path: str, local_path: str
    ) -> None:
        """
        Function to upload files/folders to workspaces.
        Args:
              workspace_id (int) : id of the where file need to uploaded
              workspace_path (str) : file path on workspace if folder does not exist it will create
              local_path (str) : uploaded file path
              
        Raises:
              InvalidParameterException: when the parameter like workspace id is invalid
              InvalidPathException: when the file to path is invalid
        
        """
  

    def download_from_workspaces(self, workspace_id: int, workspace_path: str) -> None:
        """
        Function to download files/folders from workspaces.
        Args:
              workspace_id (int) : Id of the where file need to uploaded
              workspace_path (str) : Downloaded file on workspace
        Returns:
              None
        Raises:
              InvalidPathException : for invalid path
              OperationFailedException : when downloading fails
              InvalidParameterException: when the parameter like workspace id is invalid
       
        """
      
