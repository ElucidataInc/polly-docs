class jobs:
    """
    The polly_jobs class contains functions which can be used to create, cancel and monitor polly jobs.

    Polly CLI jobs can now be initiated, managed and have a status-checked for from Polly Python.
    This lets users run jobs on the Polly cloud infrastructure by scaling computation resources as per need.
    Users can start and stop tasks and even monitor jobs.

    Args:
        token (str): token copy from polly.

    Usage:
        from polly.jobs import jobs
        
        jobs = jobs(token)
    """

    def submit_job(self, project_id: str, job_file: str) -> pd.DataFrame:
        """
        Submits  a polly cli job in the given workspace
        
        Args:
            project_id (str) :  workspace id
            job_file (str) : required configuration json file path
            
        Returns:
            on success, a dataframe with workspace id and job id
            on failure, throws errors
        """

    def cancel_job(self, project_id: str, job_id: str):
        """
        Cancel a polly job.
        
        Args:
            project_id (str) : workspace id
            job_id (str) : job id
        
        Raises:
            InvalidParameterException
        """


    def job_status(self, project_id: str, job_id="", internalCalls=False) -> dict:
        """
        Get the status of a job given the rproject id and job id.
        If no job id give, gives status for all the jobs in the
        provided workspace
        
        Args:
            project_id (str) : workspace id
            job_id (str) : job id
            internalCalls :
        
        Returns:
            dataframe with job id, job name and status sorted as per created timestamp
        """
