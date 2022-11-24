class Cohort:
    """
The Cohort class contains functions which can be used to create cohorts, add or remove samples, \
merge metadata and data-matrix of samples/datasets in a cohort and edit or delete a cohort.

Args:
	token -- token copy from polly: (str)
    
Examples:
	>>> from polly.cohort import Cohort
          
	>>> cohort = Cohort(token)
	If you are authorised then you can initialize object without token.
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
            This function is used to create a cohort.
            Args:
                  local_path(str): local path to instantiate the cohort.
                  cohort_name(str): identifier name for the cohort.
                  description(str): description about the cohort.
                  repo_key(str): Optional argument: repo_key(repo_name/repo_id) for the omixatlas from where \
datasets or samples is to be added.
                  entity_id(list): Optional argument: list of dataset_id or sample_id to be added to the \
cohort.
            Returns:
                  A confirmation message on creation of cohort.
            Errors:
                  InvalidParameterException: Empty or Invalid Parameters
                  InvalidCohortNameException: The cohort_name does not represent a valid cohort name.
                  InvalidPathException: Provided path does not represent a file or a directory.
             After making Cohort Object you can create cohort.
             
	Examples:
             1. while passing argument of repo and dataset.
              >>> cohort.create_cohort("/import/tcga_cohort","cohort_name","cohort_description",\
"repo_id",list_of_datasets)

             2. without passing argument of repo and dataset.
              >>> cohort2.create_cohort("/path","cohort_name","cohort_description")
        """

    def add_to_cohort(self, repo_key: str, dataset_id=None, sample_id=None) -> None:
        """
        This function is used to add datasets or samples to a cohort.
        Args:
              repo_key(str): repo_key(repo_name OR repo_id) for the omixatlas where datasets or samples belong.
              entity_id(list): list of dataset ID or sample ID to be added to the cohort.
        Returns:
              A confirmation message for number of datasets or samples which are added to the cohort.
        Errors:
              InvalidParameterException: Empty or Invalid Parameters.
              InvalidCohortOperationException: This operation is not valid as no cohort has been \
instantiated.
         After creating cohort we can add datasets or samples to cohort.
         Example-
        .. code::
                cohort.add_to_cohort("repo_id",list_of_dataset_ids)
        """

    def remove_from_cohort(self, dataset_id=None, sample_id=[]) -> None:
        """
        This function is used for removing datasets or samples from a cohort.
        Args:
              entity_id(list): list of dataset IDs or sample IDs which is to be removed from the cohort.
        Returns:
              A confirmation message on removal of datasets or samples from cohort.
        Errors:
              InvalidParameterException: Empty or Invalid Parameters
              InvalidCohortOperationException: This operation is not valid as no cohort has been \
instantiated.
          Example-
        .. code::
                cohort.remove_from_cohort(list_of_datasets)
        """

    def merge_data(self, data_level: str):
        """
        Function to merge metadata (dataset, sample and feature level metadata) or data-matrix of all the samples\
or datasets in the cohort.
        Args:
            | data_level(str): identifier to specify the data to be merged - "dataset", "sample", "feature" or \
"data_matrix"
        Returns:
            | A pandas dataframe containing the merged data which is ready for analysis
        """

    def is_valid(self) -> bool:
        """
        This function is used to check the validity of a cohort.
        Returns:
            |  A boolean result based on the validity of the cohort.
        Errors:
            |  InvalidPathException: Cohort path does not represent a file or a directory.
            |  InvalidCohortOperationException: This operation is not valid as no cohort has been \
instantiated.
        .. code::
                cohort2.is_valid()
        """

    def load_cohort(self, local_path: str):
        """
        Function to load an existing cohort into an object. Once loaded, the functions described in the \
documentation can be used for the object where the cohort is loaded.
        Args:
            |  local_path(str): local path of the cohort.
        Returns:
            |  A confirmation message on instantiation of the cohort.
        Errors:
            |  InvalidPathException: This path does not represent a file or a directory.
            |  InvalidCohortPathException: This path does not represent a Cohort.
        |  Example-
        .. code::
                cohort.load_cohort("/path/cohort_name.pco")
        """

    def edit_cohort(self, new_cohort_name=None, new_description=None):
        """
          This function is used to edit the cohort level metadata such as cohort name and description.
        Args:
            |  new_cohort_name(str): Optional Argument: new identifier name for the cohort.
            |  new_description(str): Optional Argument: new description about the cohort.
        Returns:
            |  A confirmation message on updation of cohort.
        Errors:
            |  InvalidCohortOperationException: This operation is not valid as no cohort has been \
instantiated.
            |  CohortEditException: No parameter specified for editing in cohort
        |  Example-
        .. code::
                cohort.edit_cohort("edited-cohort-name","edited-cohort-description")
        """

    def summarize_cohort(self):
        """
        Function to return cohort level metadata and dataframe with datasets or samples added in the cohort.
        Returns:
              A tuple with the first value as cohort metadata information (name, description and number of \
dataset(s) or sample(s) in the cohort) and the second value as dataframe containing the source, \
dataset_id or sample_id and data type available in the cohort.
        Errors:
              InvalidCohortOperationException: This operation is not valid as no cohort has been \
instantiated.
         Example-
        .. code::
                metadata, cohort_details = cohort.summarize_cohort()
         metadata will contain a object like this.
        .. code::
                {
                'cohort_name': 'cohort_name',
                'number_of_samples': 6,
                'description': 'cohort_description'
                }
         cohort detail will contain a table like that
        .. csv-table::
            :header: "", "source_omixatlas",  "datatype", "dataset_id"
            :delim: |
            0 |	tcga |	Mutation |	BRCA_Mutation_TCGA-A8-A09Q-01A-11W-A019-09
            1 |	tcga |	Mutation |	BRCA_Mutation_TCGA-AN-A0FL-01A-11W-A050-09
            2 |	tcga |	Mutation |	BRCA_Mutation_TCGA-AR-A254-01A-21D-A167-09
            3 |	tcga |	Mutation |	BRCA_Mutation_TCGA-D8-A1XO-01A-11D-A14K-09
            4 |	tcga |	Mutation |	BRCA_Mutation_TCGA-EW-A1J2-01A-21D-A13L-09
            5 |	tcga |	Mutation |	BRCA_Mutation_TCGA-LL-A50Y-01A-11D-A25Q-09
        """
        
    def delete_cohort(self) -> None:
        """
        This function is used to delete a cohort.
        Returns:
            | A confirmation message on deletion of cohort
        """

    def create_merged_gct(self, file_path: str, file_name="") -> None:
        """
          This function is used to merge all the gct files in a cohort into a single gct file.
        Args:
            |  file_path(str): the system path where the gct file is to be written.
            |  file_name(str): Optional Argument: Identifier for the merged file name, cohort name will be used by default.
        """
