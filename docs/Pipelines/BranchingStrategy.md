
Your feature branch must be forked off from master ideally. Suppose, you are creating a new pipeline named pipeline_007 pertaining to Jira ticket TKT-162, you should create a branch like this:

```
git checkout master
git checkout -b TKT-162/ftr/pipeline_007_dev
```
<br>

There are 3 main branches in the pipelines repository - `develop`, `staging` & `master`. Your branch should be first merged to `develop` followed by `staging` and then when it is production-ready, a last PR should be raised to `master`. Each branch merge or commit to 3 main branches or on `*_dev` branches triggers a deployment to specific environment(s). Check the following table to understand deployment triggers:

| Branch    | Devpolly          | Testpolly         | Polly             | Description                                                                   |
| --------- | ----------------- | ----------------- | ----------------- | ----------------------------------------------------------------------------- |  
| `*_dev`   | :material-check:  | :material-close:  | :material-check:  | Commits to branches ending with `_dev` are deployed to devpolly and on polly  |
| `develop` | :material-check:  | :material-close:  | :material-check:  | PR merges to `develop` branch is deployed to devpolly and polly               |
| `staging` | :material-close:  | :material-check:  | :material-check:  | PR merges to `staging` branch is deployed to testpolly and polly              |
| `master`  | :material-close:  | :material-close:  | :material-check:  | PR merges to `master` branch is only deployed to polly                        |


!!! note "Important Note"
    As you can see, all the 4 branch types (`*_dev`, `develop`, `staging`, `master`) are deployed to production along with other environments. This approach ensures that pipeline developers can test their pipelines without relying solely on the stability of devpolly and testpolly environments. In the production environment, each pipeline is assigned a 'stage' attribute that differentiates its maturity level. The *_dev and develop branches are designated as the 'dev' stage, the staging branch represents the 'test' stage, and the master branch is considered the 'prod' or production stage.

<br>
<br>
<br>
<br>
<br>
<br>

