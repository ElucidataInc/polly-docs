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
            repo_key (str): repo_id OR repo_name. This is a mandatory field.
            body (dict): The payload is a JSON file which should be as per the structure defined for\
            schema. Only data-admin will have the authentication to update the schema
        
        Raises:
            apiErrorException: Some Issue in Replacing the Schema for the OmixAtlas.
            paramException: Parameter Functions are not passed correctly.

        Returns:
            dict: returns the replaced schema
        
        """

    def update_schema(self, repo_key: str, body: dict) -> dict:
        """Update Schema for an existing OmixAtlas

        Args:
            repo_key (str): repo_id OR repo_name
            body (dict): The payload is a JSON file which should be as per the structure defined for\
                schema. Only data-admin will have the authentication to update the schema.

        Raises:
            apiErrorException: Some Issue in Inserting the Schema for the OmixAtlas.
            paramException: Parameter Functions are not passed correctly.

        Returns:
            dict: Dict containing the Updated Schema
        """

    def insert_schema(self, repo_key: str, body: dict) -> dict:
        """Insert the Schema for an OmixAtlas

        Args:
            repo_key (str): repo_id OR repo_name.
            body (dict):  The payload is a JSON file which should be as per the structure defined for \
                schema. Only data-admin will have the authentication to update the schema.

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
        """ Function to Get Schema For an OmixAtlas

        Args:
            repo_key (str): repo_id OR repo_name. This is a mandatory field
            schema_level (list, optional): The default value is all the table names for the repo.
            Defaults to [].Users can also a specific table name on which they want to query the schema.
            Users can table names using `SHOW TABLES IN <repo>` query.
            Also backward compatible with previous schema_level values of ['dataset', 'sample'].
            source (str, optional): is the source from where data is ingested into the OmixAtlas. 
                                    Defaults to "".
            The default value is 'all', which will fetch the schema of all sources.
            data_type (str, optional): is the datatype for which user wants to get the schema for.
            The default value is 'all', which will fetch the schema of all datatypes. Defaults to "".
            return_type (str, optional): Defaults to "dataframe". Users can also get output in dict
                                         format

        Raises:
            paramException: When Function Parameter passed are not in the right format.
            RequestException: There is some issue in fetching the Schema
            invalidApiResponseException: The Data returned is not in appropriate format

        Returns:
            dict, dataframe: It will contain the schema for specific table names as dict or dataframe depending on return_type
        """
