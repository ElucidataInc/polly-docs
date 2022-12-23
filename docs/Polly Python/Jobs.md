
::: polly.Jobs
    options:
      show_source: false

## Examples

``` python
import os
from polly.auth import Polly
from polly.jobs import jobs
token = os.environ['POLLY_REFRESH_TOKEN']
Polly.auth(token)
job = jobs()
```

``` python
job_map = {
    "cpu":"1",
    "memory":"1Gi",
    "image": "ubuntu",
    "tag": "latest",
    "name": "This is a job for testing bust mode",
    "command": ["/bin/bash","-c", "TERM=xterm free -h; echo '\nnumber of vCPU';nproc;sleep 30"]
  }
```
``` python
import json
json_object = json.dumps(job_map, indent=4)
with open ("/import/test_job.json", "w") as input:
    input.write(json_object)
```
### submit_job()

``` python
#job submission
workspace_id = "14832"
job_file = "/import/test_job.json"
job.submit_job(workspace_id,job_file)
```

``` 
   Workspace ID                            Job ID
0         14832  1f445778c0034b44a799a0e72d48ab94
```

### job_status()

``` python
#getting job status immediately after running 
job.job_status(workspace_id,"1f445778c0034b44a799a0e72d48ab94")
```


``` 
                             Job ID                             Job Name  \
0  1f445778c0034b44a799a0e72d48ab94  This is a job for testing bust mode   

  Job State  
0   RUNNING  
```

``` python
#post job completion
job.job_status(workspace_id,"1f445778c0034b44a799a0e72d48ab94")
```
``` 
                             Job ID                             Job Name  \
0  1f445778c0034b44a799a0e72d48ab94  This is a job for testing bust mode   

  Job State  
0   Success  
```

``` python
#getting status of all jobs
job.job_status(workspace_id)
```


``` 
                              Job ID                             Job Name  \
0   cc1ad0d41eb64cc68dd12c4ddbefad9a  This is a job for testing bust mode   
1   9d9895b453fe42fcb13bb7902e5de457  This is a job for testing bust mode   
2   d9e7c71c6a154a2bb018458c07db0733  This is a job for testing bust mode   
3   6d9617e623444ec2888ded77c1affa2a  This is a job for testing bust mode   
4   2e1a2d8ac504446ea081d40dcc7ca1f1  This is a job for testing bust mode   
5   af9934a82dfe4207a85d7c59e62d9f7b  This is a job for testing bust mode   
6   ff21afe991b641a689d93df801c2329e  This is a job for testing bust mode   
7   b87ffa99f7ba469c8a34ddd17d647290  This is a job for testing bust mode   
8   2f235fe366304ef580a2893059a122ee  This is a job for testing bust mode   
9   9d27a3e3e1b242799f614889b81f28b1  This is a job for testing bust mode   
10  8298ddd6895242cc80c4f490b8f11d17  This is a job for testing bust mode   
11  fcb749ed5af941aebe21b87c1f963bb1  This is a job for testing bust mode   
12  b5e7c5317b7040e2b4da86334e9e12d4  This is a job for testing bust mode   
13  a4822e0d7c624f679d93ef84a412abb4  This is a job for testing bust mode   
14  90ee553e98b54246909a826cc5569727  This is a job for testing bust mode   
15  01dd2ca6e2884668872f087fdae7f8ca  This is a job for testing bust mode   
16  67bbed35a210421c8988ab30e975037e  This is a job for testing bust mode   
17  01c931994e1f4a71aac17d895735ee4f  This is a job for testing bust mode   
18  496277121f0045e18e6e1a9a0c4005db  This is a job for testing bust mode   
19  5a45874082e14e57ad6989e38911e6f1  This is a job for testing bust mode   
20  ed8664bc71f94a87bba8c0ee99c5690a  This is a job for testing bust mode   
21  1ae5cc5a32934cd2815c04a879b82c8e  This is a job for testing bust mode   
22  1660474c0f54452b8bd802403d01c9ac  This is a job for testing bust mode   
23  edd1ee1cf1d946e5a1d5f03d38f1a933  This is a job for testing bust mode   
24  d259389229304cd585487bf114008f2c  This is a job for testing bust mode   
25  8cc98c884e6e422baf117508bc47ef2b  This is a job for testing bust mode   
26  972d75d64d4c40889a540eb792f214c4  This is a job for testing bust mode   
27  11d6c4e70230418aa047ab2cc3266917  This is a job for testing bust mode   
28  c4a70af3ae7544658d2aad1f9a4364d2  This is a job for testing bust mode   
29  813918d6f5644441929012d8b29954eb  This is a job for testing bust mode   
30  87fd6e4c44bc4294922d8351492c4316  This is a job for testing bust mode   
31  aa6f922b3529455a976e72e51f220fbf  This is a job for testing bust mode   
32  b1282ca0454a4c3da3e4d304ac602876  This is a job for testing bust mode   
33  b083e9f7ad204071a4ad3b71005e6f99  This is a job for testing bust mode   
34  5fd2065c5ef241a89ea8d0cfdfb1e7a0  This is a job for testing bust mode   
35  73a0650927b542d396b8199fbf08d85b  This is a job for testing bust mode   
36  a4c7a7b08f82421c90158d826b8c1f0a  This is a job for testing bust mode   
37  7695077a47fb41358ca14c291b932b8e  This is a job for testing bust mode   
38  6483eda8bfef4d52b22cbcf9318efe78  This is a job for testing bust mode   
39  3449e1f8925d41b98a48a00069253c27  This is a job for testing bust mode   
40  ec00c7b67eaf4816af5cb6f3e28ab842  This is a job for testing bust mode   
41  9dc9dd97871d49f7a733862b9aec7ba0  This is a job for testing bust mode   
42  1b3c7b3916a84120b103a4da18a53c47  This is a job for testing bust mode   
43  79060530f30441a29c87e1d8447734ec  This is a job for testing bust mode   
44  d1f61d879c6e4d06a2bb4fe518da61eb  This is a job for testing bust mode   
45  ff1d4a55155b475aadcaa28158a04607  This is a job for testing bust mode   
46  8e9982ea324f4abfaf5f62566fff2947  This is a job for testing bust mode   
47  103a51688f7b4132938e644123486c56  This is a job for testing bust mode   
48  a94f01d9c0084d9a8c2a3fcebda2e001  This is a job for testing bust mode   
49  535ff212c8fa4d3484f5561afc67ef45  This is a job for testing bust mode   
50  5ce824808aa9465d8c023621fe993d78  This is a job for testing bust mode   
51  963cc0883c0a4bf092037eb9b1568477  This is a job for testing bust mode   
52  02c5e93a5a74486dbe5efe9288fdcd8c  This is a job for testing bust mode   
53  1c761cb0064c4710a74a50f524150d5f  This is a job for testing bust mode   
54  b71d62d020df441aac7986e4e5d32186  This is a job for testing bust mode   
55  a3bbebd0cdef491e811b0c2b3b8852cc  This is a job for testing bust mode   
56  3b220d0a324444ae805a00ce26b36d00  This is a job for testing bust mode   
57  5842dc069db64296b3212c3b7935e80d  This is a job for testing bust mode   
58  1f445778c0034b44a799a0e72d48ab94  This is a job for testing bust mode   

    Job State  
0     Success  
1     Success  
2     Success  
3     Success  
4     Success  
5     Success  
6   CANCELLED  
7   CANCELLED  
8     Success  
9     Success  
10    Success  
11    Success  
12    Success  
13  CANCELLED  
14  CANCELLED  
15  CANCELLED  
16    Success  
17  CANCELLED  
18  CANCELLED  
19  CANCELLED  
20  CANCELLED  
21  CANCELLED  
22  CANCELLED  
23  CANCELLED  
24  CANCELLED  
25  CANCELLED  
26    Success  
27    Success  
28  CANCELLED  
29  CANCELLED  
30    Success  
31    Success  
32    Success  
33    Success  
34    Success  
35    Success  
36    Success  
37    Success  
38    Success  
39    Success  
40    Success  
41    Success  
42  CANCELLED  
43    Success  
44    Success  
45    Success  
46    Success  
47    Success  
48    Success  
49    Success  
50    Success  
51    Success  
52    Success  
53    PENDING  
54    Success  
55    Success  
56    Success  
57    Success  
58    Success  
```

``` python
#workspace id can be int or string
workspace_id = 14832
job.job_status(workspace_id, "1f445778c0034b44a799a0e72d48ab94")
```
``` python
# invalid workspace id 
workspace_id = 148324345
job.job_status(workspace_id, "1f445778c0034b44a799a0e72d48ab94")
```

    Not able to get the status of the Job(s)
    Unexpected exception formatting exception. Falling back to standard exception

``` 
                             Job ID                             Job Name  \
0  1f445778c0034b44a799a0e72d48ab94  This is a job for testing bust mode   

  Job State  
0   Success  
```

``` python
# invalid job id provided while status
workspace_id = 14832
job.job_status("14832", "1f445778c0034b44a799a0e72d4adac")
```
``` 
Not able to get the status of the Job(s)
Unexpected exception formatting exception. Falling back to standard exception
``` 

### cancel_job()

``` python
#cancelling a job post success 
# "{'errors': [{'status': '400', 'code': 'bad_req', 'title': 'Bad Request', 'detail': 'Cannot cancel a job in Success state'}]}" -> this will be removed 
job.cancel_job(workspace_id, "1f445778c0034b44a799a0e72d48ab94")
```
``` 
{'errors': [{'status': '400', 'code': 'bad_req', 'title': 'Bad Request', 'detail': 'Cannot cancel a job in Success state'}]}
Failed to cancel the job.: Cannot cancel a job in Success state
``` 

``` python
#cancelling jobs before success
job.submit_job(workspace_id,job_file)
job.cancel_job(workspace_id, "a51187d93aea402b8f40789404398bfb")
```

``` 
Cancelled job ID a51187d93aea402b8f40789404398bfb successfully!
``` 

``` python
#cancelling already cancelled job
# {'errors': [{'status': '400', 'code': 'bad_req', 'title': 'Bad Request', 'detail': 'Cannot cancel a job in CANCELLED state'}]} -> this will be removed 
job.cancel_job(workspace_id, "a51187d93aea402b8f40789404398bfb")
```

``` 
{'errors': [{'status': '400', 'code': 'bad_req', 'title': 'Bad Request', 'detail': 'Cannot cancel a job in CANCELLED state'}]}
Failed to cancel the job.: Cannot cancel a job in CANCELLED state
``` 

``` python
# cancelling job with invalid job id 
job.cancel_job(workspace_id, "a51187d93aea402b8f40789404398bmn")
```
``` 
---------------------------------------------------------------------------
    InvalidJobFunctionParameterException      Traceback (most recent call last)
    Input In [22], in <cell line: 2>()
          1 # cancelling job with invalid job id 
    ----> 2 job.cancel_job(workspace_id, "a51187d93aea402b8f40789404398bmn")
    
    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:120, in jobs.cancel_job(self, project_id, job_id)
        117     raise InvalidParameterException("Missing/invalid datatype for job id")
        119 if not self._is_valid_job_id(project_id, job_id):
    --> 120     raise InvalidJobFunctionParameterException("Job id")
        122 self.base_url = f"https://v2.api.{self.session.env}.elucidata.io"
        123 self.jobUrl = f"{self.base_url}/projects/{project_id}/jobs/{job_id}"
    
InvalidJobFunctionParameterException: The specified Job id could not be found. Inspect and try again.
``` 

``` python
# cancelling job with invalid workspace id 
job.cancel_job(12345, "a51187d93aea402b8f40789404398bfb")
```
``` 
Not able to get the status of the Job(s)
Unexpected exception formatting exception. Falling back to standard exception
``` 
### job_logs()

```py
# logs are not generated yet 
job.job_logs("14832","a8e49e1d545c493f9cddf60f90baac2d")
```

    Error: Not able to find the logs. It seems to be not generated yet

```py
# parameter checks -> incorrect job id 
job.job_logs("14832","a8e49e1d545c493f9cddf60fzzzzzz")
```

    Error: Not able to get the log of the job



    ---------------------------------------------------------------------------

    InvalidJobFunctionParameterException      Traceback (most recent call last)

    Input In [24], in <cell line: 2>()
          1 # parameter checks -> incorrect job id 
    ----> 2 job.job_logs("14832","a8e49e1d545c493f9cddf60fzzzzzz")


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:274, in jobs.job_logs(self, project_id, job_id, mode)
        272 except Exception as err:
        273     print("Error: Not able to get the log of the job")
    --> 274     raise err


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:237, in jobs.job_logs(self, project_id, job_id, mode)
        235 parameter_check_dict["job_id"] = job_id
        236 parameter_check_dict["mode"] = mode
    --> 237 self._parameter_check_for_jobs(parameter_check_dict)
        238 if isinstance(project_id, int):
        239     project_id = str(project_id)


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:597, in jobs._parameter_check_for_jobs(self, paramter_dict)
        595 elif keys == "job_id":
        596     job_id = values
    --> 597     self._check_job_id(job_id, project_id)
        598 elif keys == "mode":
        599     mode = values


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:609, in jobs._check_job_id(self, job_id, project_id)
        607     raise InvalidParameterException("Missing/invalid datatype for job id")
        608 if not self._is_valid_job_id(project_id, job_id):
    --> 609     raise InvalidJobFunctionParameterException("Job id")


    InvalidJobFunctionParameterException: The specified Job id could not be found. Inspect and try again.



```py
# parameter checks -> incorrect mode
job.job_logs("14832","a8e49e1d545c493f9cddf60f90baac2d","kjb,yjv")
```

    Error: Not able to get the log of the job



    ---------------------------------------------------------------------------

    InvalidParameterException                 Traceback (most recent call last)

    Input In [26], in <cell line: 2>()
          1 # parameter checks -> incorrect mode
    ----> 2 job.job_logs("14832","a8e49e1d545c493f9cddf60f90baac2d","kjb,yjv")


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:274, in jobs.job_logs(self, project_id, job_id, mode)
        272 except Exception as err:
        273     print("Error: Not able to get the log of the job")
    --> 274     raise err


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:237, in jobs.job_logs(self, project_id, job_id, mode)
        235 parameter_check_dict["job_id"] = job_id
        236 parameter_check_dict["mode"] = mode
    --> 237 self._parameter_check_for_jobs(parameter_check_dict)
        238 if isinstance(project_id, int):
        239     project_id = str(project_id)


    File /usr/local/lib/python3.10/site-packages/polly/jobs.py:601, in jobs._parameter_check_for_jobs(self, paramter_dict)
        599 mode = values
        600 if mode not in DATA_LOG_MODES:
    --> 601     raise InvalidParameterException(
        602         f"valid mode values are : {DATA_LOG_MODES}"
        603     )


    InvalidParameterException: Empty or Invalid Parameters = valid mode values are : ['latest', 'all'].



```py
# parameter checks -> no job id provided
job.job_logs("14832")
```


    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    Input In [27], in <cell line: 2>()
          1 # parameter checks -> incorrect mode
    ----> 2 job.job_logs("14832")


    TypeError: jobs.job_logs() missing 1 required positional argument: 'job_id'
