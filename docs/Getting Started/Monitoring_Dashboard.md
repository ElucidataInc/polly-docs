## What is a Monitoring Dashboard?

Polly's Monitoring Dashboard will allow users to track and monitor all ingestions and runs made to the Omix Atlases mapped to your Org. This tool allows you to keep track of all your ingestion runs, view logs in the event of a failure, and categorize your runs.

All of your ingestions, whether they are running, completed or errored, are displayed when you first access the dashboard.

### Ingestion runs:

Each run displayed on the dashboard will include information such as,

1. Run name - OA to which ingestion was performed, timestamp of when the run began, and number of datasets ingested, in the order respectively

For eg:



2. Status - whether run is pending, in progress, completed or errored.
3. Priority - The urgency with which an ingestion run is introduced depends on its priority.
4. User name - this essentially refers to the users who started the ingestion run.
5. Duration - Displays the amount of time that has passed since the run was initiated.
6. Actions - The action button will lead you to the details page

### Filters:

To select and categorize your runs, utilize these five filters available.

1. It has two options - My Runs/All Runs
<ul>
 <li>**My runs** - is a filter that allows all ingestion runs in the user's pipeline to be displayed.
 <li>**All runs** - view all ingestion runs performed by the organization's employees on the OAs to which the user has access
    - The level of access will determine the visibility of the OAs ingestions
</ul></li>
2. The Priority of the Ingestion Runs are categorized as **High, Medium,** and **Low** in the second filter. The urgency with which an ingestion run will be consumed depends on its priority. It will be easier to prepare for potential delays if you are aware of the ingestion run's priority.

3. Third, Filter has bifurcation for **In progress Runs and All runs that are completed or stopped, or in progress**.
4. The fourth filter will classify the results according to status:

- **Completed** - This is the optimal status, which means that all tasks in the ingestion run were successfully completed.
- **In-Progress** - indicates that one or more jobs have begun running and that the run is literally in progress. Some tasks may have also been completed.
- **Errored** - This state denotes that at least one or more ingestion run jobs failed before they could complete; the dataset has otherwise been processed, and nothing is ongoing.
- **Pending** - This status indicates that the run is still awaiting ingestion to begin (only been registered in-queue)

5. The last filter is to use the calendar to list the runs according to when they were triggered.

### View the Details page:

The top section shows the run name, the number of datasets being added to the OA, who started the run, and the amount of time that has passed ( since the run started )

Then, at the bottom half of the page, there are two sections: one is titled "Dataset Progress," and the other is "Logs."

- **Dataset Progress tab** - Users can get a list of all the datasets in a run along with their progress status under the Progress tab. The run's overall progress can also be seen. It consists of three parts: the dataset ID, the accompanying dataset's status, and the amount of progress made.

### To check the logs:

Click on the log tab

There are three items under the logs tab: a search option, a button to download all logs, and a part where logs are presented.

- **Search option** - Users can look for a certain dataset or sample using the search feature in the logs.
- **Download all logs button** - This will be enabled only when the run is completed. The entire detailed logs will be accessible for download for improved accessibility and usability.once the run is finished so that users can inspect, explore, and determine the source of any errors.
- **Main section of logs** - The logs will be updated here once each minute. The user will have access to a thorough record that will allow them to track and keep tabs on the progress of their processes.
