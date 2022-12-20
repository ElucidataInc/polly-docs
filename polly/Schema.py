class Schema:
    """
The OmixAtlas class contains functions to enable users to interact with the schema of a particular OmixAtlas
Functions for visualizing, updating and inserting schema is released.
Updating and inseting schema is allowed for users who have data-admin credentials only.

Args:
    token (str): token copy from polly.
    
Usage:
    from polly.OmixAtlas import OmixAtlas
    omixatlas = OmixAtlas(token)
    """

    def replace_schema(self, repo_key: str, body: dict) -> dict:
        """Replace the same for an Omixatlas. The function will replace the existing source and datatype
        dictionary if new source and datatype dictionaries are passed

        Args:
            repo_key (str): repo_id OR repo_name of the OmixAtlas. Users can get this by running get_all_omixatlas function.
            body (dict): The payload should be a JSON file for a specific table as per the structure defined for schema. 
            Only data-admin at organisation level will have the authentication to replace the schema.
        
        Raises:
            apiErrorException: Some Issue in Replacing the Schema for the OmixAtlas.
            paramException: Parameter Functions are not passed correctly.

        Returns:
            dict: returns the replaced schema
        
        """

    def update_schema(self, repo_key: str, body: dict) -> dict:
        """Update Schema for an existing OmixAtlas

        Args:
            repo_key (str): repo_id OR repo_name of the OmixAtlas. Users can get this by running get_all_omixatlas function.
            body (dict): The payload should be a JSON file for a specific table as per the structure defined for schema. 
            Only data-admin at organisation level will have the authentication to update the schema.

        Raises:
            apiErrorException: Some Issue in Inserting the Schema for the OmixAtlas.
            paramException: Parameter Functions are not passed correctly.

        Returns:
            dict: Dict containing the updated Schema
        """

    def insert_schema(self, repo_key: str, body: dict) -> dict:
        """ 
        This function is used to insert the Schema in a newly created OmixAtlas.

        Args:
            repo_key (str): repo_id OR repo_name of the OmixAtlas. Users can get this by running get_all_omixatlas function.
            body (dict):  The payload should be a JSON file for a specific table as per the structure defined for schema. 
            Only data-admin at organisation level will have the authentication to insert the schema.

        Raises:
            apiErrorException: Some Issue in Inserting the Schema for the OmixAtlas.

        Returns:
            dict: Dict containing the New Schema
        """
    def get_schema(
        self,
        repo_key: str,
        schema_level=[],
        source="",
        data_type="",
        return_type="dataframe",
    ) -> dict:
        """ 
        Function to get the Schema of all the tables in an OmixAtlas. Any user can get the schema of an OmixAtlas.

        Args:
             repo_key (str): repo_id OR repo_name. This is a mandatory field.
             schema_level (list, optional): Table name for which users want to get the schema. \
             Users can get the table names by querying `SHOW TABLES IN <repo_name>` using query_metadata function.\
             The default value is all the table names for the repo.
             source (str, optional): Source for which user wants to fetch the schema. 
             The default value is all the sources in the schema.
             data_type (str, optional): Datatype for which user wants to fetch the schema. 
             The default value is all the datatypes in the schema.
             return_type (str, optional): For users who intend to query should use "dataframe" output. 
             For users, who want to perform schema management, they should get the output in "dict" format. 
             Dataframe format doesn't give the complete schema, it only shows the information \
             which aids users for writing queryies. Default value is "dataframe". 

        Raises:
            paramException: When Function Parameter passed are not in the right format.
            RequestException: There is some issue in fetching the Schema
            invalidApiResponseException: The Data returned is not in appropriate format

        Returns:
            dict, dataframe: It will contain the schema for specific table names as dict or dataframe depending on return_type
        """
