It can’t be stressed enough how important Data is. And so protecting it against any kind of theft is of utmost priority. We follow the standard norms and best practices to ensure customer data is safe with us. Steps we take to ensure data safety and integrity.

*  Private Databases 
Databases are one of our primary data storage places. We have taken the necessary steps to make them accessible only by a trusted entity. This is achieved by hosting our databases in a Virtual private cloud which essentially makes them impossible to be accessed directly via the internet.

*  Encryption at rest 
Data uploaded by users on our platform, for analysis or storage purposes is encrypted using server-side encryption. At the time of computation, all disk volumes used for ephemeral storage are encrypted using an industry-standard AES 256 encryption.

*  Encryption in transit 
Data remains in transit phase from the time the user uploads the data to the time it finally gets stored on our platform. For this duration, the data is encrypted & hence protected using TLS 2.1 protocol.


## Preventing Unauthorized access

Data breaches due to compromised & weak passwords is one of the major reasons why data gets landed into the wrong hands. Data privacy is also a major concern. We don’t want to allow anyone to get access to our private data unless explicitly allowed. To make sure data accessing and operations on that data are performed by someone who owns it or is allowed to do so, we follow the following measures

* Principle of least privilege

* Role-based access control

* MFA

* Strong Password protection via AWS Cognito


## Infrastructure security

* **Network Isolation:** All AWS compute resources, responsible for processing and storing of user data, are in a separate Virtual Private Cloud(VPC). VPC’s are isolated networks within the AWS cloud

* **Private Dockers**: All Docker Images pushed to Polly’s Docker Domain are private by default and are only accessible by system admins 



## AWS Security promises

AWS promises to take care of the infrastructure and its security. Below is the list of the compliances and services that follow those compliances

* **SOC**(1,2,3), **PCI**, **HIPAA**

* [AWS VPC](https://aws.amazon.com/vpc/ "https://aws.amazon.com/vpc/")

* [AWS S3](https://aws.amazon.com/s3/ "https://aws.amazon.com/s3/")

* [AWS Cloudwatch](https://aws.amazon.com/cloudwatch/ "https://aws.amazon.com/cloudwatch/") & [Cloudwatch logs](https://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/WhatIsCloudWatchLogs.html "https://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/WhatIsCloudWatchLogs.html")

* [AWS Cognito](https://aws.amazon.com/cognito/ "https://aws.amazon.com/cognito/")

* [RDS](https://aws.amazon.com/rds/ "https://aws.amazon.com/rds/")

* [DynamoDB](https://aws.amazon.com/dynamodb/ "https://aws.amazon.com/dynamodb/") 



For a complete list of compliances that AWS services adhere to, please head over [here](https://aws.amazon.com/compliance/services-in-scope/ "https://aws.amazon.com/compliance/services-in-scope/") 


## Security whitepaper

For a more detailed report, please read our [security whitepaper](https://ss-usa.s3.amazonaws.com/c/308473941/media/103395fa8e7b4f0c0568268917274591/Security%20and%20Compliance%20in%20Polly%20-%20Accelerating%20drug%20discovery%20secur.pdf "https://ss-usa.s3.amazonaws.com/c/308473941/media/103395fa8e7b4f0c0568268917274591/Security%20and%20Compliance%20in%20Polly%20-%20Accelerating%20drug%20discovery%20secur.pdf")
