class Download:
    """
The OmixAtlas class contains functions to download gct, h5ad, and vcf files, and the metadata of any dataset.
    """
    def download_data(self, repo_id, _id: str):
        """
        To download any dataset, the following function can be used to get the signed URL of the dataset.
        The data can be downloaded by clicking on this URL. 
        NOTE: This signed URL expires after 60 minutes from when it is generated.
         
        The repo_name OR repo_id of an OmixAtlas can be identified by calling the get_all_omixatlas() function.
        The dataset_id can be obtained by querying the metadata at the dataset level using query_metadata().
        This data can be parsed into a data frame for better accessibility using the code under the examples section.

        Args:
              repo_id (str/int): repo_id for the omixatlas
              payload (dict): The payload is a JSON file which should be as per the structure defined for schema.
              Only data-admin will have the authentication to update the schema.
           
        Raises:
              apiErrorException: Params are either empty or its datatype is not correct or see detail.
        """

    def download_metadata(self, repo_key: str, dataset_id: str, file_path: str) -> None:
        """
        This function is used to download the dataset level metadata into a json file.
        
        Args:
            repo_key (str): repo_key(repo_name/repo_id) of the repository which is linked.
            dataset_id(str): dataset_id of the dataset to be linked.
            file_path(str): the system path where the json file is to be written.
       
        """
