
# Polly disaster recovery - Mar 2021

All the data being ingested and/or generated inside Polly is being stored inside AWS. Types of data include input/output files, logs, DB transactional data & static data. No data as of now is being stored outside of the AWS account. Although AWS provides very high durability and reliability for all the data stores, there is still a risk associated with some circumstances like accidental deletion, natural hazard & unwanted access to the AWS account.

<br/>

## Datastores

There are majorly 5 data stores in which the data is being stored currently inside AWS:

-   [**RDS**](https://aws.amazon.com/rds/ "https://aws.amazon.com/rds/") **(Postgres)**
    
-   [**DynamoDB**](https://aws.amazon.com/dynamodb/ "https://aws.amazon.com/dynamodb/")
    
-   [**Elasticsearch**](https://aws.amazon.com/elasticsearch-service/ "https://aws.amazon.com/elasticsearch-service")
    
-   [**EFS**](https://aws.amazon.com/efs/ "https://aws.amazon.com/efs/")
    
-   [**S3 buckets**](https://aws.amazon.com/s3/ "https://aws.amazon.com/s3/")

<br/>    

|Datastore|Purpose|Backup|Frequency|
|--|--|--|--
|RDS (Postgres)|**For Tabular data** <br/> To store relational transactional data, including but not limited to users & organizations' information, workspaces data/metadata etc.|Uses [AWS backup service](https://aws.amazon.com/backup/ "https://aws.amazon.com/backup/") (managed service for automated backups)|Daily|
|DynamoDB|**For JSON based data** <br/> It works as extension to Postgres for the core working of the platform.|Uses [AWS backup service](https://aws.amazon.com/backup/ "https://aws.amazon.com/backup/") (managed service for automated backups)|Daily|
|Elasticsearch|**For Analysis Ready data**|Before ingesting data into Elasticsearch, it is being loaded to an Amazon S3 bucket(s). In case Elasticsearch data is lost, we can re-index the data from S3 buckets.|Daily|
|EFS|**For application state storing** <br/> Applications based on shiny architecture uses EFS for storing application state at a given point in time.|Uses [AWS backup service](https://aws.amazon.com/backup/ "https://aws.amazon.com/backup/") (managed service for automated backups)|Daily|
|S3 buckets|**For File storage** <br/><br/> S3 contains: <br/> **Static content:** Media, Configurations etc <br/>**Input/Output files:** Files ingested into the platform and/or generated output files <br/> **Logs** Application logs|Custom cron job exectution <br/><br/>The job copies the files every day into a separate account.|Daily|

<br/>
<br/>

## Key take-aways

-   Backups are done at least once a day
    
-   Following datastores are backed up - [**RDS**](https://aws.amazon.com/rds/ "https://aws.amazon.com/rds/")**,** [**DynamoDB**](https://aws.amazon.com/dynamodb/ "https://aws.amazon.com/dynamodb/")**,** [**Elasticsearch**](https://aws.amazon.com/elasticsearch-service/ "https://aws.amazon.com/elasticsearch-service/")**,** [**EFS**](https://aws.amazon.com/efs/ "https://aws.amazon.com/efs/")**,** [**S3 buckets**](https://aws.amazon.com/s3/ "https://aws.amazon.com/s3/")**.**
    
-   **In case of a disaster, a maximum of one day of data will be lost.**