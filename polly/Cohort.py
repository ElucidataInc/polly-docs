class Cohort:
    """
The Cohort class contains functions which can be used to create cohorts, add or remove samples, \
merge metadata and data-matrix of samples/datasets in a cohort and edit or delete a cohort.

Args:
	token (str): token copy from polly
    
Usage:
	from polly.cohort import Cohort
          
	cohort = Cohort(token)
    """

    def create_cohort(
        self,
        local_path: str,
        cohort_name: str,
        description: str,
        repo_key=None,
        dataset_id=None,
        sample_id=None,
    ) -> None:
        """
            This function is used to create a cohort. After making Cohort Object you can create cohort.
	    
            Args:
                  local_path (str): local path to instantiate the cohort.
                  cohort_name (str): identifier name for the cohort.
                  description (str): description about the cohort.
                  repo_key (str/int): repo_key(repo_name/repo_id) for the omixatlas from where datasets or samples is to be added.
                  entity_id (list): list of dataset_id or sample_id to be added to the cohort.
		  
            Returns:
                  A confirmation message on creation of cohort.
		  
            Raises:
                  InvalidParameterException: Empty or Invalid Parameters
                  InvalidCohortNameException: The cohort_name does not represent a valid cohort name.
                  InvalidPathException: Provided path does not represent a file or a directory.
        """

    def add_to_cohort(self, repo_key: str, dataset_id=None, sample_id=None) -> None:
        """
        This function is used to add datasets or samples to a cohort.
	
        Args:
              repo_key (str/int): repo_key(repo_name OR repo_id) for the omixatlas where datasets or samples belong.
              entity_id (list): list of dataset ID or sample ID to be added to the cohort.
	      
        Returns:
              A confirmation message for number of datasets or samples which are added to the cohort.
        
	Raises:
              InvalidParameterException: Empty or Invalid Parameters.
              InvalidCohortOperationException: This operation is not valid as no cohort has been instantiated.
        """

    def remove_from_cohort(self, dataset_id=None, sample_id=[]) -> None:
        """
        This function is used for removing datasets or samples from a cohort.
	
        Args:
              entity_id (list): list of dataset IDs or sample IDs which is to be removed from the cohort.
	      
        Returns:
              A confirmation message on removal of datasets or samples from cohort.
	      
        Raises:
              InvalidParameterException: Empty or Invalid Parameters
              InvalidCohortOperationException: This operation is not valid as no cohort has been instantiated.
        """

    def merge_data(self, data_level: str):
        """
        Function to merge metadata (dataset, sample and feature level metadata) or data-matrix of all the samples or datasets in the cohort.
        
	Args:
		data_level (str): identifier to specify the data to be merged - "dataset", "sample", "feature" or "data_matrix"
        
	Returns:
		A pandas dataframe containing the merged data which is ready for analysis
        """

    def is_valid(self) -> bool:
        """
        This function is used to check the validity of a cohort.
        
	Returns:
		A boolean result based on the validity of the cohort.
        
	Raises:
		InvalidPathException: Cohort path does not represent a file or a directory.
		InvalidCohortOperationException: This operation is not valid as no cohort has been instantiated.
        """

    def load_cohort(self, local_path: str):
        """
        Function to load an existing cohort into an object. 
	Once loaded, the functions described in the documentation can be used for the object where the cohort is loaded.
	
        Args:
		local_path (str): local path of the cohort.
		
        Returns:
		A confirmation message on instantiation of the cohort.
        Raises:
		InvalidPathException: This path does not represent a file or a directory.
		InvalidCohortPathException: This path does not represent a Cohort.
        """

    def edit_cohort(self, new_cohort_name=None, new_description=None):
        """
	This function is used to edit the cohort level metadata such as cohort name and description.
        Args:
		new_cohort_name (str): Optional Argument: new identifier name for the cohort.
		new_description (str): Optional Argument: new description about the cohort.
        
	Returns:
		A confirmation message on updation of cohort.
        
	Raises:
		InvalidCohortOperationException: This operation is not valid as no cohort has been instantiated.
		CohortEditException: No parameter specified for editing in cohort
        """

    def summarize_cohort(self):
        """
        Function to return cohort level metadata and dataframe with datasets or samples added in the cohort.
	
        Returns:
              A tuple with the first value as cohort metadata information (name, description and number of dataset(s) or sample(s) in the cohort) and the second value as dataframe containing the source, dataset_id or sample_id and data type available in the cohort.
        
	Raises:
              InvalidCohortOperationException: This operation is not valid as no cohort has been instantiated.
        """
        
    def delete_cohort(self) -> None:
        """
        This function is used to delete a cohort.
	
        Returns:
		A confirmation message on deletion of cohort
        """

    def create_merged_gct(self, file_path: str, file_name="") -> None:
        """
        This function is used to merge all the gct files in a cohort into a single gct file.
	
        Args:
		file_path (str): the system path where the gct file is to be written.
		file_name (str): Identifier for the merged file name, cohort name will be used by default.
        """
