class Basic:
    """
   OmixAtlas class enables users to interact with functional properties of the omixatlas such as \
   create and update an Omixatlas, get summary of it's contents, add, insert, update the schema, \
   add, update or delete datasets, query metadata, download data, save data to workspace etc.\

Args:
    token (str): token copy from polly.
    
Usage:
    from polly.OmixAtlas import OmixAtlas
    omixatlas = OmixAtlas(token)
    """

    def get_all_omixatlas(self, query_api_version="v2", count_by_source=True, count_by_data_type=True):
        
        """
        This function will return the summary of all the Omixatlas on Polly which the user has access to.
        Args:
              None
        Returns:
            It will return a list of JSON objects. (See Examples)
        """

    def omixatlas_summary(
        self,
        repo_key: str,
        query_api_version="v2",
        count_by_source=True,
        count_by_data_type=True):
         
        """
        This function will return you a object that contain summary of a given Omixatlas.
        Args:
            repo_key (str/int): repo_key(repo_name/repo_id)
        Returns:
            It will return a JSON object. (see examples)
        """

    def create(
        self,
        display_name: str,
        description: str,
        repo_name="",
        image_url="",
        components=[]) -> pd.DataFrame:
        
        """
        This function is used to create a new omixatlas.
        Args:
             display_name (str): display name of the omixatlas.
             description (str): description of the omixatlas.
             repo_name (str): repo_name for that Omixatlas
             image_url (str): Url of the icon for omixatlas. Optional Parameter.
             initials (str): Initials shown in the icon of omixatlas. Optional Parameter.
             explorer_enabled (bool): Default True. Optional Parameter.
             studio_presets (list): Optional Paramter.
             components (list): Optional Parameter.
        Returns:
             Dataframe after creation of omixatlas.
        Raises:
              ValueError: Repository creation response is in Incorrect format.

        """

    def update(
        self,
        repo_key: str,
        display_name="",
        description="",
        image_url="",
        components=[],) -> pd.DataFrame:
        """
        This function is used to update an omixatlas.
        Args: 
             repo_key (str/int): repo_key(repo_name/repo_id) for that Omixatlas
             display_name (str): display name of the omixatlas. Optional Parameter
             description (str): description of the omixatlas. Optional Parameter
             image_url (str): Url of the icon for omixatlas. Optional Parameter
             components (list): List of components to be added. Optional Parameter
        """
        
