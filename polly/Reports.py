class Reports:
    """
    The OmixAtlas class contains functions which will enable users to generate reports,
    link reports to a dataset in OmixAtlas and fetch reports linked with dataset in an OmixAtlas.

    Args:
        token (str): token copy from polly.

    Usage:
        from polly.OmixAtlas import OmixAtlas
        
        omixatlas = OmixAtlas(token)
    """
    def generate_report(self, repo_key: str, dataset_id: str, workspace_id: int, workspace_path="") -> None:
       """
       This function is used to generate a metadata distribution report for a dataset belonging to the geo repository.
       The generated report is then saved to the workspace provided.
       This is a MVP release to minimise time taken by users to determine relevance of a dataset in their research,
       we’re enabling auto-generation of reports.
       These reports will contain dataset and sample level metadata along with some graphs showing how the samples are distributed.
       It will help them figure out which cohorts could be of interest in a given dataset.
       This template can be modified and we’re open to user’s feedback.
       This report is available on the user’s local path as well as an option to upload the report to workspaces is given.
 
       Args:
           repo_key (str): repo_key(repo_name/repo_id) for which the report is to be generated
           Dataset_id (str): dataset_id for which the report is to be generated
           Workspace_id (str): workspace_id to where the report is to be uploaded
           workspace_path (str): workspace_path to which the report is to be uploaded
 
       Raises:
           InvalidParameterException: Empty or Invalid Parameters
           UnsupportedRepositoryException: repo other than GEO/ Unsupported repo
           UnauthorizedException: unauthorized to perform this task
       """

    def link_report(
       self,
       repo_key: str,
       dataset_id: str,
       workspace_id: int,
       workspace_path: str,
       access_key: str) -> None:
       """
       This function is used to link a file (html or pdf) present in a workspace with the specified dataset in OmixAtlas.       
       On success it displays the access key URL and a success message.
       Org admins can now link a report present in a workspace with the specified datasets in an OmixAtlas.
       Once a report is linked to a dataset in OmixAtlas, it can be fetched both from front-end and polly-python.
 
       Args:
           repo_key (str): repo_name/repo_id of the repository to be linked
           dataset_id (str): dataset_id of the dataset to be linked
           workspace_id (str): workspace_id for the file which is to be linked
           workspace_path (str): workspace_path for the file which is to be linked
           access_key (str): access_key (private or public) depending upon the link access type to be generated. 
           If public, then anyone with a Polly account with the link will be able to see the report. 
           If private, then only the individuals who have access to the workspace where reports is stored will be able to see them.
 
       Raises:
           InvalidParameterException: invalid parameters
        
       Returns:
           None
       """
    def fetch_linked_reports(self, repo_key: str, dataset_id: str) -> pd.DataFrame:
       """
       Fetch linked reports for a dataset_id in an omixatlas. 
 
       Args:
           repo_key (str): repo_key(repo_name/repo_id) of the repository for which to fetch the report
           Dataset_id (str): dataset_id of the dataset which to fetch the reports.
 
       Raises:
           InvalidParameterException : invalid parameters
           RequestException : api request exception
           UnauthorizedException : unauthorized to perform this action
 
       Returns:
           A Dataframe with the details of the linked reports - who added the report, when it was added and the link.

       """ 
    def delete_linked_report(self, repo_key: str, dataset_id: str, report_id: str) -> None:
       """
       Delete the link of the report in workspaces with the specified dataset in OmixAtlas. On success displays a success message.
 
       Arguments:
           repo_key (str): repo_key(repo_name/repo_id) of the repository which is linked.
           Dataset_id (str): dataset_id of the dataset to be unlinked
           Report_id (str): report id associated with the report in workspaces that is to be deleted. This id can be found when invoking the fetch_linked_report() function.
 
       Raises:
           InvalidParameterException: Invalid parameter passed
      
       Returns:
           None
       """       

