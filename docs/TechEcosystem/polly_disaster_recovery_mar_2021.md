All the data being ingested and/or generated inside Polly is being stored inside AWS. Types of data include input/output files, logs, DB transactional data & static data. No data as of now is being stored outside of the AWS account. Although AWS provides very high durability and reliability for all the data stores, there is still a risk associated with some circumstances like accidental deletion, natural hazard & unwanted access to the AWS account.

## Data-stores

There are majorly 5 data stores in which the data is being stored currently inside AWS:

*   [**RDS (Postgres)**](https://aws.amazon.com/rds/ "https://aws.amazon.com/rds/") 
    
*   [**DynamoDB**](https://aws.amazon.com/dynamodb/ "https://aws.amazon.com/dynamodb/")
    
*   [**Elasticsearch**](https://aws.amazon.com/elasticsearch-service/ "https://aws.amazon.com/elasticsearch-service")
    
*   [**EFS**](https://aws.amazon.com/efs/ "https://aws.amazon.com/efs/")
    
*   [**S3 buckets**](https://aws.amazon.com/s3/ "https://aws.amazon.com/s3/")

<br/>    

|Datastore|Purpose|Backup|Frequency|
|--|--|--|--
|RDS (Postgres)|**For Tabular data** <br/> To store relational transactional data, including but not limited to users & organizations' information, workspaces data/metadata etc.|Uses [AWS backup service](https://aws.amazon.com/backup/ "https://aws.amazon.com/backup/") (managed service for automated backups)|Daily|
|DynamoDB|**For JSON based data** <br/> It works as extension to Postgres for the core working of the platform.|Uses [AWS backup service](https://aws.amazon.com/backup/ "https://aws.amazon.com/backup/") (managed service for automated backups)|Daily|
|Elasticsearch|**For Analysis Ready data**|Before ingesting data into Elasticsearch, it is being loaded to an Amazon S3 bucket(s). In case Elasticsearch data is lost, we can re-index the data from S3 buckets.|Daily|
|EFS|**For application state storing** <br/> Applications based on shiny architecture uses EFS for storing application state at a given point in time.|Uses [AWS backup service](https://aws.amazon.com/backup/ "https://aws.amazon.com/backup/") (managed service for automated backups)|Daily|
|S3 buckets|**For File storage** <br/><br/> S3 contains: <br/> **Static content:** Media, Configurations etc <br/>**Input/Output files:** Files ingested into the platform and/or generated output files <br/> **Logs** Application logs|Custom cron job exectution <br/><br/>The job copies the files every day into a separate account.|Daily|

<br/>

## Measures for Disaster Recovery and Management

To ensure we don't lose customers’ data, we have the following contingency plan in place.

**Automated backups**

Daily Backup of persistent storage is performed for both SQL and NoSQL databases and EFS Filesystems. Recovery has been given key importance here, to make sure the downtime can be reduced to a few hours in case of a serious unwanted event. We use AWS backup services for Dynamo, RDS, EFS, EBS, and S3 storage services. In instances where AWS backup services are not available / feasible (eg. Opensearch) we have a custom in-house support process in place. We have a 120 days point in time restoration which allows us to restore any desired data within this time period. The data is encrypted at rest and has a cross regional backup. 

**Versioning of files** 

Versioning of files is enabled for S3 to recover accidental file deletions upon a users' request.



## FAQs

**1. How does Elucidata ensure the security of Polly's infrastructure and data?**

Along with holding SOC 2 and AWS FTR approvals, Elucidata enforces strict security protocols, conducts regular internal audits, and ensures all staff are trained in security awareness. This comprehensive approach helps secure the infrastructure and data of Polly.

**2. How does Polly ensure protection against known vulnerabilities in system components and software?**

Polly extensively uses serverless architecture components like AWS API Gateway, Lambda, S3, CloudFront, etc., managed by AWS, where Elucidata does not need to explicitly manage patching. For services like EC2 instances, patches are applied quarterly. Polly's patching cadence and vulnerabilities are monitored in real-time using securityscorecard.com, with identified issues typically remediated within two weeks.

**3. What measures are in place for data protection law compliance?**

Our platform prioritizes your privacy and security. While we don’t store any Personal Identifiable Information (PII) or Protected Health Information (PHI), we ensure streamlined and secure login processes by retaining only essential customer details such as names and email addresses.

 **4. Are all connections to Polly secured?**

Yes, all connections to Polly are secured via SSL/TLS protocols, ensuring integrity and security against unauthorized access.

**5.How is data backup and encryption managed?**

Backups are encrypted and stored in a separate AWS region. They are performed daily and retained for 60 days.

**6. How does Polly ensure personal information is stored securely?**

All Personal Information is encrypted using strong cryptographic methods. AWS Key Management Service with AES-256 in GCM mode is employed for encryption.
