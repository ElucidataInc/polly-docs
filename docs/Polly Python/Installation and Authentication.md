## Install Polly Python using pip

<pre><code>pip install polly-python</code></pre>

## Import from libraries

The following libraries need to be imported over the development environment to access the data.

<pre><code>from polly.auth import Polly
from polly.omixatlas import OmixAtlas
from polly.workspaces import Workspaces
from polly.cohort import Cohort</code></pre>

## Authentication
Authentication of the account is required to be able to access the capabilities of the Polly Python library.

### Copying authentication token or key
To get this token, follow the steps below:

1. Go to [Polly](https://polly.elucidata.io)<br>

2. Click the **User Options** icon from the left-most panel<br>

3. Click on **Authentication** on the panel that appears<br>

4. Click on **Copy** to copy the authentication token or key<br>

### Using the token
The following code is required to add the authentication function in the Polly Python library

<pre><code>AUTH_TOKEN = "[authentication_token_copied]"
Polly.auth(AUTH_TOKEN)</code></pre>

### Using the key
Go to authentication keys and input the name, description to generate a key. This can be copied and used as shown below:-

<pre><code>AUTH_KEY = "[authentication_key_copied]"
Polly.auth(AUTH_KEY)</code></pre>
