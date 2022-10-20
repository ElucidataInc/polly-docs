## Polly Workflow

### Why Nextflow:

-   Significant reduction in execution time with parallel processing of pipeline/sub-processes within a pipeline with support of polly's computationally resources

-   Offer users with interactive reports, timelines, DAGs, and trace documents—which we don't have currently but which will help with monitoring and boost the pipeline efficiency through better planning.

-   Reduced effort and more portability as dockerization of all the codes is not needed

### Nextflow on Polly:

- All Polly multi-processing jobs that have a diverse machine needs and demand a lot of computing time can be converted into Nextflow pipeline, for an effective run.
- As the code uses modularized processes, it is more adaptable since individual modules may be plugged into and used in other pipelines, increasing its reproducibility and flexibility.
- For computationally intensive analyses to process enormous numbers of data and metadata, resource optimization is a major bonus, which in turn helps save cost.

---

### Comparison of Polly CLI Jobs with Nextflow on Polly :

| Features                                                                                                         | Polly Jobs | Nextflow Jobs |
|------------------------------------------------------------------------------------------------------------------|:----------:|:-------------:|
| Ready to use Parallel processing                                                                                 |     No     |      Yes      |
| Machine configuration customization for each step/process/module in a workflow/pipeline                          |     No     |      Yes      |
| Dividing pipeline into individual modules which helps in easy customization and reproducibility                  |     No     |      Yes      |
| Need for dockerization                                                                                           |     Yes    |       No      |
| Informative and interactive report documents with timelines, DAGs and resource usage stats for pipeline executed |     No     |      Yes      |

