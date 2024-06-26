Polly is a cloud based Software as a Service (SaaS) solution. Below details out various layers of Polly.

**Data Harmonization & Ingestion**

  - **Polly Pipelines:** Workflow orchestration system(supports Nextflow)
    
  - **Polly Curation:** Manual and automated data curation

**Data Storage & Management**

  - **Polly Atlas:** For storing structured data in a data model
    
  - **Polly Workspaces:** For storing unstructured data

**Data Consumption**

  - **Polly Apps:** Custom apps running on a highly scalable infra (Hail, IGV)
    
  - **Polly Notebooks:** Data exploration/analysis using Jupyter notebooks



## Infrastructure

Polly’s infrastructure is powered entirely by managed cloud providers - AWS. Backed by AWS massive infrastructure, Polly is fast, secure, highly available & scales seamlessly.

It uses managed cloud services like [AWS S3,](https://aws.amazon.com/s3/ "https://aws.amazon.com/s3/") [DynamoDB](https://aws.amazon.com/dynamodb/ "https://aws.amazon.com/dynamodb/"), [RDS](https://aws.amazon.com/rds/ "https://aws.amazon.com/rds/") (running PostgreSQL), [AWS ElasticSearch Service](https://aws.amazon.com/elasticsearch-service/ "https://aws.amazon.com/elasticsearch-service/"), [AWS Lambda](https://aws.amazon.com/lambda/ "https://aws.amazon.com/lambda/"), [Amazon API Gateway](https://aws.amazon.com/api-gateway/ "https://aws.amazon.com/api-gateway/"), [Amazon Cloudfront](https://aws.amazon.com/cloudfront/ "https://aws.amazon.com/cloudfront/"), [Amazon Cognito](https://aws.amazon.com/cognito/ "https://aws.amazon.com/cognito/"). Polly’s computational infrastructure is powered by [Kubernetes](https://kubernetes.io/ "https://kubernetes.io/") used to orchestrate [Docker](https://www.docker.com/ "https://www.docker.com/") containers on both AWS. It is highly scalable and compatible with all cloud providers.


### Micro-services & Serverless architecture

Various parts of Polly are built and deployed as micro services. This allows for better fault isolation and resilience. A focus on managed serverless components like API Gateway, Lambda functions, S3, DynamoDB, S3 etc also allows Polly to scale seamlessly under high load.

Polly’s micro-services architecture pattern allows the development teams to choose different technologies that are best suited to those parts. Polly has been built ground up using Python3.6, Angular11, NodeJS, R & R-Shiny.

### Third Party Integrations

File uploads from Polly front-end is powered by [Filestack](https://www.filestack.com/ "https://www.filestack.com/"). [Sentry](https://sentry.io/welcome/ "https://sentry.io/welcome/") integration allows Polly developers to monitor errors and fix them faster. [Mixpanel](https://mixpanel.com/ "https://mixpanel.com/") is used for gathering usage analytics which help product owners refine and improve platform capabilities and user journeys.
