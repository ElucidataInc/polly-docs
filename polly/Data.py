class Data:
    """
The OmixAtlas class contains functions to add, update, or delete datasets in an omixatlas, or save datasets from polly notebooks to workspace.
    """
    
    def dataset_metadata_template(
        self, repo_key, source="all", data_type="all"
    ) -> dict:
        """
        This function is used to create a template for dataset level metadata.

        In order to ingest the dataset level metadata appropriately in the OmixAtlas,
        the user needs to ensure the json files contains the keys as per the dataset level schema of the OmixAtlas.

        Args:
            repo_key (str/int): repo_key(repo_name/repo_id) for that Omixatlas
            source (all) : Source/Sources present in the schema. By default all
            data_type (all) : Datatype/Datatypes present in the schema. By default all
            
        """

    def save_to_workspace(
        self, repo_key: str, dataset_id: str, workspace_id: int, workspace_path: str
    ) -> json:
        """Function to save data from OmixAtlas to workspaces

        Args:
            repo_id (str): repo_id of the Omixatlas
            dataset_id (str): dataset id that needs to be saved
            workspace_id (int): workspace id in which the dataset needs to be saved
            workspace_path (str): path where the workspace resides

        Returns:
            json: Info about workspace where data is saved and of which Omixatlas
        """

    def add_datasets(
        self,
        repo_key: int,
        source_folder_path: dict,
        destination_folder_path="",
        priority="low",
    ) -> pd.DataFrame:
        """This function is used to add data to an omixatlas

        Args:
            repo_id (str/int): repo_id for that Omixatlas
            source_folder_path (dict): source folder path from data and metadata files are fetched.
            destination_folder_path (str, optional): Destination folder structure in s3. Defaults to "".
            priority (str, optional): Priority at which this data has to be inserted. Defaults to "low".

        Raises:
            paramError: If Params are not passed in the desired format or value not valid.
            RequestException: If there is issue in data ingestion.

        Returns:
            pd.DataFrame: DataFrame showing Upload Status of Files
        """

    def update_datasets(
        self,
        repo_key: int,
        source_folder_path: dict,
        destination_folder_path="",
        priority="low",
    ) -> pd.DataFrame:
        """This function is used to update data/metadata to an omixatlas

        Args:
            repo_id (str/int): repo_id for that Omixatlas
            source_folder_path (dict): source folder path from data and metadata files are fetched.
            destination_folder_path (str, optional): Destination folder structure in s3. Defaults to "".
            priority (str, optional): Priority at which this data has to be inserted. Defaults to "low".

        Raises:
            paramError: If Params are not passed in the desired format or value not valid.
            RequestException: If there is issue in data ingestion.

        Returns:
            pd.DataFrame: DataFrame showing Upload Status of Files
        """

    def delete_datasets(self, repo_key: int, dataset_ids: list):
        """This function is used to delete data from an omixatlas

        Args:
            repo_id (str/int): repo_id for that Omixatlas
            dataset_ids (list): dataset_ids that users want to delete

        Raises:
            paramError: If Params are not passed in the desired format or value not valid.
            RequestException: If there is issue in data ingestion.
        """
