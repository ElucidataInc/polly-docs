It can’t be stressed enough how important Data is. And so protecting it against any kind of theft is of utmost priority. We follow the standard norms and best practices to ensure customer data is safe with us. Steps we take to ensure data safety and integrity.

## Data Privacy & Protection

Polly’s infrastructure is hosted on AWS. All AWS compute resources are in a separate Virtual Private Cloud (VPC). VPC’s are logically isolated sections of the AWS cloud which we can define as per our requirements of a virtual cloud. We have complete control over your virtual networking environment, including the selection of your own IP address range, creation of subnets, and configuration of route tables and network gateways.

**- The principle of least privilege** : All Resources across our microservices follow the Principle of Least Privilege(POLP). This states that any resource will have only the bare minimum amount of access needed for it to function adequately. The same principle also applies to users who have access to the production systems at Elucidata barring system admins who are responsible to keep these systems in check. 

**- Encryption at rest** : All data at rest in distributed file systems or in ephemeral storage volumes during processing is encrypted using AES 256 encryption. Also, all other implicit storage volumes are encrypted using AWS encryption keys 

**- Encryption in transit** : Protects your data if communications are intercepted while data moves between your site and the cloud provider or between two services. This protection is achieved by encrypting the data before transmission; authenticating the endpoints, and decrypting and verifying the data on arrival. We have encrypted all data transfers with TLS. 

**- Database security** : Relational Databases have been hidden behind a firewall and only accessible from within the same Virtual private cloud or by system admin from an already provisioned bastion host. Only necessary inbound traffic from other microservices to the database is allowed and SSH ports are closed. 

**- Private docker domain** : All Docker Images pushed to Polly’s Docker Domain are private by default. No other user, either internal or external has access to these Docker images barring system admins. Upon a mail request, we can delete the data within 60 days.


## Measures for authorization, authentication, and accountability

Polly has been designed by embracing Zero Trust Policy as one of its core security principles. Several safety measures have been put in place even for cases where an attacker might have access to compromised credentials.

**- Multi-Factor Authentication (MFA):** All Administrators on the production Account for AWS have Multi-Factor Authentication(MFA) enabled to mitigate the possibility of a break-in.
  
**- Access logs:** Upon the users' request and use case, we can provide you with access to AWS logs for both resources and users who had access to the production environment.

**- Encrypted passwords:** We do not store users' passwords in our database. They are salted and stored within AWS Cognito wherein no member of Elucidata has the ability to view them. Encryption of secrets like STS tokens(temporary authorisation tokens) is done at rest.

**- Role-based access control:** To enforce the Principle of Least Privilege Access(POLP) in a more effective way, resources, and operations are managed by a robust Role-based access control system (RBAC). The system allows an admin in an organization to create teams, map Polly components to these teams and assign data access levels depending on user roles.

**- Health monitoring:** We use prometheus to monitor our kubernetes clusters, resource usage (RAM, CPU and disk), pods, nodes and kubernetes services. We have also configured alerts for unexpected events on the PollyHealth Dashboard, which is used to learn about the availability and operation status of Polly’s services. AWS is responsible for protecting Polly's infrastructure under their shared Responsibility Model.

**- Risk assessment:** We do Third Party VAPT (Vulnerability Assessment and Penetration Testing) audits twice in a financial year. Vulnerability Assessments are done every quarter. We also perform third party risk assessment annually to monitor and mitigate risks across the organization.

## Infrastructure security

* **Network Isolation:** All AWS compute resources, responsible for processing and storing of user data, are in a separate Virtual Private Cloud(VPC). VPC’s are isolated networks within the AWS cloud

* **Private Dockers**: All Docker Images pushed to Polly’s Docker Domain are private by default and are only accessible by system admins 


## Security whitepaper

For a more detailed report, please read our [security whitepaper](../img/Architecture/Security_whitepaper.pdf)



## FAQs

**1. Are third parties given access to any proprietary information from Polly?**
  
  No, proprietary information from Polly is not shared with third-party vendors.

**2. Are Polly's security policies reviewed regularly?**

  Yes, Polly's Information Security Policy is reviewed once a year to ensure it remains robust and up-to-date.

**3. Does Elucidata have a dedicated security team? How many members, and do you work with external vendors?**

  Elucidata has a dedicated internal security team, including a VP of Engineering, Director of Engineering, Sr. Manager of Operations, IT Manager, and Cybersecurity Specialist. We also engage with external vendors like Akitra Inc., Astra Security, Prescient Assurance, SentinelOne, and SpringVerify for various security functions.

**4. Does Polly monitor and log access to network resources?**

Yes, Polly uses AWS Config and CloudTrail to create an audit trail for its infrastructure and network.

**5. What is the password policy for the customer application infrastructure within Polly?**

  Passwords must be at least 8 characters with no upper length limit, containing upper and lowercase alphabetic characters, one numeric 
  character, and one special character. Users can change passwords at their will from within the Polly GUI, and password rotation is not mandated at predefined intervals. Passwords are encrypted and stored in AWS Cognito, ensuring compliance with several industry standards, and are encrypted in transit using TLS 1.3.

**6. How does Polly ensure personal information is stored securely?**

  All Personal Information is encrypted using strong cryptographic methods. AWS Key Management Service with AES-256 in GCM mode is employed for encryption.


## AWS Security promises

AWS promises to take care of the infrastructure and its security. Below is the list of the compliances and services that follow those compliances

* **SOC(1,2,3)**, **PCI**, **HIPAA**

* [AWS VPC](https://aws.amazon.com/vpc/ "https://aws.amazon.com/vpc/")

* [AWS S3](https://aws.amazon.com/s3/ "https://aws.amazon.com/s3/")

* [AWS Cloudwatch](https://aws.amazon.com/cloudwatch/ "https://aws.amazon.com/cloudwatch/") & [Cloudwatch logs](https://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/WhatIsCloudWatchLogs.html "https://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/WhatIsCloudWatchLogs.html")

* [AWS Cognito](https://aws.amazon.com/cognito/ "https://aws.amazon.com/cognito/")

* [RDS](https://aws.amazon.com/rds/ "https://aws.amazon.com/rds/")

* [DynamoDB](https://aws.amazon.com/dynamodb/ "https://aws.amazon.com/dynamodb/") 



For a complete list of compliances that AWS services adhere to, please head over [here](https://aws.amazon.com/compliance/services-in-scope/ "https://aws.amazon.com/compliance/services-in-scope/") 



