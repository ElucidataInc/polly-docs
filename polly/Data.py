class Data:
    """
The OmixAtlas class contains functions to add, update, or delete datasets in an omixatlas, or save datasets from polly notebooks to workspace.
    """
    
    def dataset_metadata_template(
        self, repo_key, source="all", data_type="all"
    ) -> dict:
        """
        This function is used to fetch the template of dataset level metadata in a given OmixAtlas. 
        In order to ingest the dataset level metadata appropriately in the OmixAtlas, 
        the user needs to ensure the metadata json files contains the keys as per the dataset level schema of the OmixAtlas.


        In order to ingest the dataset level metadata appropriately in the OmixAtlas,
        the user needs to ensure the json files contains the keys as per the dataset level schema of the OmixAtlas.

        Args:
            repo_key (str/int): repo_key(repo_name/repo_id) for that Omixatlas
            source (all, optional) : Source/Sources present in the schema. Default value is "all"
            data_type (all, optional) : Datatype/Datatypes present in the schema. Default value is "all"
            
        """

    def save_to_workspace(
        self, repo_id: str, dataset_id: str, workspace_id: int, workspace_path: str
    ) -> json:
        """
        Function to download a dataset from OmixAtlas and save it to Workspaces.

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
        repo_id: int,
        source_folder_path: dict,
        destination_folder_path="",
        priority="low",
    ) -> pd.DataFrame:
        """
        This function is used to add a new data into an OmixAtlas. Once user runs this function successfully, 
        they should be able to see the ingestion status on the data ingestion monitoring dashboard after ~15 mins. 

        Args:
            repo_id (str/int): repo_id for that Omixatlas
            source_folder_path (dict): source folder paths from data and metadata files are fetched. In this dictionary, there should be two \
            keys called "data" and "metadata" with value consisting of folders where data and metadata is stored respectively.  
            destination_folder_path (str, optional): Destination folder structure in s3. 
            Users should use this only when they want to manage the folder structure in the backend. 
            It is advised to not not give any value for this, by default the data goes in root folder.
            priority (str, optional): Priority at which this data has to be ingested into the OmixAtlas. 
            The default value is "low". Acceptable values are "medium" and "high".

        Raises:
            paramError: If Params are not passed in the desired format or value not valid.
            RequestException: If there is issue in data ingestion.

        Returns:
            pd.DataFrame: DataFrame showing Upload Status of Files
        """

    def update_datasets(
        self,
        repo_id: int,
        source_folder_path: dict,
        destination_folder_path="",
        priority="low",
    ) -> pd.DataFrame:
        """
        This function is used to update an existing data into an OmixAtlas. Once user runs this function successfully, 
        they should be able to see the ingestion status on the data ingestion monitoring dashboard after ~15 mins. 

        Args:
            repo_id (str/int): repo_id for that Omixatlas
            source_folder_path (dict): source folder paths from data and metadata files are fetched. In this dictionary, there should be two \
            keys called "data" and "metadata" with value consisting of folders where data and metadata is stored respectively.
            destination_folder_path (str, optional): Destination folder structure in s3. Users should use this only when \
            they want to manage the folder structure in the backend. It is advised to not not give any value for this, by default the data goes in \
            root folder.
            priority (str, optional): Priority at which this data has to be ingested into the OmixAtlas. The default value is "low". \
            Acceptable values are "medium" and "high".

        Raises:
            paramError: If Params are not passed in the desired format or value not valid.
            RequestException: If there is issue in data ingestion.

        Returns:
            pd.DataFrame: DataFrame showing Upload Status of Files
        """

    def delete_datasets(self, repo_id: int, dataset_ids: list):
        """
        This function is used to delete datasets from an OmixAtlas. Once user runs this function successfully, 
        they should be able to see the ingestion status on the data ingestion monitoring dashboard after ~15 mins. 

        Args:
            repo_id (str/int): repo_id for that Omixatlas
            dataset_ids (list): list of dataset_ids that users want to delete.

        Raises:
            paramError: If Params are not passed in the desired format or value not valid.
            RequestException: If there is issue in data ingestion.
        """
