Cohort class of polly-python enables users to create cohorts, add/remove datasets or samples from them, merge the dataset, sample, feature and data-matrix level metadata across all samples of the cohort, delete a cohort etc.
This feature is enabled in both types of OmixAtlas:- 

where 1 dataset has 1 sample such as TCGA, Depmap, LINCS, cBioportal, CPTAC, Immport and GDC.

where 1 dataset has multiple samples such as GEO, Metabolomics OmixAtlas etc.

In order to enable cohorting on a given OmixAtlas, please contact polly.support@elucidata.io

### create_cohort()
Cohort creation is enabled in the local environment - be it in the polly notebook environment or user's local. The minimum requirement for a cohort is to have a cohort.meta file inside the cohort that defines the .pco format. The cohort.meta file is encrypted in base64 format for keeping the metadata consistent and protected.
```
cohort.create_cohort(local_path=”<path>”,cohort_name=”name”,description=”description”, 
                    repo_key=”repo_key” (optional), dataset_id=list (optional))
```
### add_to_cohort()
This function allows users to add datasets to the cohort.
```
cohort.add_to_cohort(repo_key=”<repo_id or repo_name>”,dataset_id=[“dataset_id1”,…])
```
### remove_from_cohort()
This function removes the samples from a cohort. 
```
cohort.remove_from_cohort(dataset_id=[“dataset_id1”,…]))
```
### summarize_cohort()
It returns a tuple with the first value as cohort metadata information (name, description and number of dataset(s) or sample(s) in the cohort) and the second value as dataframe containing the source, dataset_id or sample_id  and data type available in the cohort.
```
cohort.summarize_cohort(dataset_id=[“dataset_id1”,…]))
```
### load_cohort()
This function loads an already existing cohort into a newly instantiated object for working on the cohort.
```
cohort.load_cohort(local_path=”path to cohort”)
```
### edit_cohort()
This feature is used for renaming cohort_name and/or cohort description from cohort level metadata.
```
cohort.edit_cohort(new_cohort_name=”new_name”,new_cohort_description=”new description”)
```
### merge_data()
Function to merge the dataset level metadata from all the GCT files in a cohort. Returns a pandas Dataframe containing the merged data for analysis.
```
cohort.merge_data("dataset")
```
Function to merge the sample level metadata from all the GCT files in a cohort. Returns a pandas Dataframe containing the merged data for analysis.
```
cohort.merge_data("sample")
```
Function to merge the feature level metadata from all the GCT files in a cohort. Returns a pandas Dataframe containing the merged data for analysis.
```
cohort.merge_data("feature")
```
Function to merge the data-matrix level metadata from all the GCT files in a cohort. Returns a pandas Dataframe containing the merged data for analysis.
```
cohort.merge_data("data_matrix")
```
### delete_cohort()
This function deletes an existing cohort.
```
cohort.delete_cohort()
```
### is_valid()
This function is for validating a cohort. This functions returns a boolean result depending on the validity of the cohort. 
```
cohort.is_valid()
```

### create_merged_gct()

This function is used to create a merged gct file from a cohort

```
from polly.cohort import Cohort
cohort = Cohort()

cohort.load_cohort(“cohort_path”)

cohort.create_merged_gct(file_path,file_name: optional)
```

Argument description:-

file_path (str): path where the merged gct file should be saved

file_name (str): (optional) file name of the merged gct. By default, the cohort name will be used

