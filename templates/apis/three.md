# <svg t="1702051717164" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="5796" width="40" height="40"><path d="M501.07 557.26m-252.29 0a252.29 252.29 0 1 0 504.58 0 252.29 252.29 0 1 0-504.58 0Z" fill="#D2DAFF" p-id="5797"></path><path d="M914.5 185.45c-2.5-24.22-7.91-39.32-17.53-48.94S872.26 121.48 848 119c-14.86-1.54-33.57-2-55.61-1.3-34.24 1-65.63 4.39-70.73 5-83.03 5.23-163.66 34.21-233.26 83.75A457.08 457.08 0 0 0 375 322c-9.52-0.48-61-1.85-118.81 19.32-58.82 21.56-135.48 71.12-166.62 184.8l46.82 21.12c21.86-33.35 47.42-52.71 76-57.55 44.93-7.61 85.69 22.24 87.89 23.88l0.81-1 0.51 2.62L518.36 732l1.67 0.32 3.81 4.88c0.31 0.4 31.48 41.83 23.94 87.52-4.73 28.68-24.12 54.36-57.61 76.33l21.12 46.81c90.11-24.68 154-83.83 184.79-171.07 20-56.6 19.67-107.49 19.18-120.81A457.84 457.84 0 0 0 827 545.08c49.55-69.57 78.52-150.23 83.81-233.29 0.57-5.1 3.92-36.49 5-70.73 0.67-22.06 0.19-40.75-1.31-55.61zM335.39 394.31l-27.08 63.36c-23.34-12.22-61.45-26.41-104.1-19.38q-4.8 0.79-9.51 1.91 31.59-32.18 77.66-49.48a274.83 274.83 0 0 1 72-16q-4.73 9.7-8.97 19.59z m311.66 365.16Q629.93 808 597.56 841q0.87-3.88 1.52-7.85c7.51-45.52-9.17-85.88-21.82-108.61l61.91-26.46q11.77-5 23.22-10.77a294.93 294.93 0 0 1-15.34 72.16zM863.67 244.4c-1.22 32.49-4.53 61.9-4.56 62.19l-0.12 1.32c-9.16 147.29-105.72 284.88-240.25 342.37l-76 32.5-192-192 32.5-76.05c57.49-134.53 195.08-231.08 342.37-240.25l1.32-0.12c0.29 0 29.7-3.34 62.19-4.56 47.33-1.78 65.1 2.14 70.63 4 1.78 5.5 5.7 23.27 3.92 70.6z" fill="#27187F" p-id="5798"></path><path d="M183.23 880.41l-30.16-30.16a615.34 615.34 0 0 1 18.31-71.72c17.14-53.63 38.14-93.13 62.42-117.41 29.13-29.12 56.52-43.12 83.72-42.85 22.91 0.25 44.39 10.72 65.68 32s31.76 42.77 32 65.68c0.31 27.2-13.71 54.59-42.85 83.73C348.08 824 308.58 845 255 862.1a613.37 613.37 0 0 1-71.77 18.31z m133.46-210.15c-12.11 0-28.05 9.54-46.13 27.62-26.74 26.75-46.23 81.53-57.27 122.31 40.78-11 95.56-30.52 122.31-57.27 18.21-18.21 27.76-34.25 27.62-46.38-0.06-5.52-1.86-14.57-16.78-29.49s-24-16.73-29.49-16.79zM620.92 499A86.43 86.43 0 1 1 682 473.67 85.89 85.89 0 0 1 620.92 499z m0-120.85a34.42 34.42 0 1 0 24.35 10.08 34.21 34.21 0 0 0-24.35-10.1z" fill="#27187F" p-id="5799"></path></svg> Sample Applications of the ADMETlab 3.0

---

We have tried to simplify the user's operation by encapsulating the computational modules as much as possible. In short,
the application scenario of this platform focuses on a single application of ADMET property prediction of compounds, so
the main steps can be abbreviated as: inputting the molecules, cleaning the molecules (optional), modeling the
computation, and outputting the results. This section lists some code snippets that we believe are commonly used by
users. If you have any questions, please feel free to contact us via the Contact page.

## 1. Wash molecules

<span style="color: #2878b5">Example code:</span>
<ul class="nav nav-tabs" id="myTab" role="tablist">
  <li class="nav-item" role="presentation">
    <button class="nav-link active" id="home-tab" data-bs-toggle="tab" data-bs-target="#python" type="button" role="tab" aria-controls="home" aria-selected="true"><i class="iconfont icon-python"></i> Python</button>
  </li>
  <li class="nav-item" role="presentation">
    <button class="nav-link" id="profile-tab" data-bs-toggle="tab" data-bs-target="#shell" type="button" role="tab" aria-controls="profile" aria-selected="false">Shell</button>
  </li>
  <li class="nav-item" role="presentation">
    <button class="nav-link" id="contact-tab" data-bs-toggle="tab" data-bs-target="#r" type="button" role="tab" aria-controls="contact" aria-selected="false">R</button>
  </li>
</ul>
<div class="tab-content" id="myTabContent">
  <div class="tab-pane fade show active" id="python" role="tabpanel" aria-labelledby="home-tab">

    import requests

    baseUrl = 'http://121.40.210.46:8097'

    if __name__ == '__main__':
        api = '/api/washmol'
        url = baseUrl + api
        param = {
            'SMILES': ["molecule", "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"],
        }
        response = requests.post(url, json=param)
        if response.status_code == 200:
            data = response.json()['data']
            print(data)

  </div>
  <div class="tab-pane fade" id="shell" role="tabpanel" aria-labelledby="profile-tab">

    #!/bin/bash

    url="http://121.40.210.46:8097/api/washmol"
    smiles="CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"

    result=$(curl -i -X POST -H "'Content-type':'application/json'" -d '{"SMILES":"'$smiles'"}' $url)
    echo $result

  </div>
  <div class="tab-pane fade" id="r" role="tabpanel" aria-labelledby="contact-tab">

    library(httr)
    library("jsonlite")

    url <- "http://121.40.210.46:8097/api/washmol"
    body <- list(
    "SMILES" = "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2",
    )

    response <- POST(url, body=body, encode="json")
    content(response)

  </div>
</div>

## 2. Predicting small data set

<span style="color: #2878b5">Example code:</span>

<ul class="nav nav-tabs" id="myTab2" role="tablist">
  <li class="nav-item" role="presentation">
    <button class="nav-link active" id="home-tab" data-bs-toggle="tab" data-bs-target="#python2" type="button" role="tab" aria-controls="home" aria-selected="true"><i class="iconfont icon-python"></i> Python</button>
  </li>
  <li class="nav-item" role="presentation">
    <button class="nav-link" id="profile-tab" data-bs-toggle="tab" data-bs-target="#input" type="button" role="tab" aria-controls="profile" aria-selected="false">Shell</button>
  </li>
</ul>
<div class="tab-content" id="myTabContent">
  <div class="tab-pane fade show active" id="python2" role="tabpanel" aria-labelledby="home-tab">

    import requests

    baseUrl = 'http://121.40.210.46:8097'

    if __name__ == '__main__':
        api = '/api/admet'
        url = baseUrl + api
        param = {
            'SMILES': ["molecule", "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"],
        }
        response = requests.post(url, json=param)
        if response.status_code == 200:
            data = response.json()['data']
            print(data)

  </div>
  <div class="tab-pane fade" id="input" role="tabpanel" aria-labelledby="profile-tab">

    import json

    import requests
    import pandas as pd

    baseUrl = 'http://121.40.210.46:8097'

    def transform(data):
        resultList = []
        for mol in data['data']:
            if not mol['data']:
                # Invalid SMILES
                tmp = {'smiles': mol['smiles']}
            else:
                tmp = dict({'smiles': mol['smiles']})
                for _, admet in mol['data'].items():
                    for endpoint in admet:
                        # endpoint is a dict
                        tmp[endpoint['name']] = endpoint['value']
            resultList.append(tmp)
        return pd.DataFrame(resultList).fillna('Invalid SMILES')

    if __name__ == '__main__':
        api = '/api/admet'
        url = baseUrl + api
        param = {
            'SMILES': []
        }
        data = pd.read_csv('example.csv')  # Read CSV file
        smiles_list = data['SMILES'].tolist()  # Getting data from the SMILES column
        param['SMILES'] = smiles_list  # Assign the SMILES list to the SMILES key of param

        response = requests.post(url, json=param)

        if response.status_code == 200:  # If access is successful
            data = response.json()['data']
            # transform to csv file
            result = transform(data)
            result.to_csv('result.csv', index=False)

  </div>
</div>

<div class="d-flex justify-content-evenly align-items-center">
    <a href="../static/apis/files/input.csv" download="input.csv">
        <button type="button" class="btn btn-secondary"><i class="bi bi-cloud-arrow-down-fill"></i> Input file</button>
    </a>
    <a href="../static/apis/files/output.csv" download="output.csv">
        <button type="button" class="btn btn-secondary"><i class="bi bi-cloud-arrow-down-fill"></i> Output file</button>
    </a>
</div>

## 3. Predicting massive amounts of molecules

For users who need to compute a large number of molecules, we recommend splitting the list of molecules into multiple
subtasks for iterative computation (1000 molecules). Submit the task and get the results in the form of code, and then
stitch the results of each subtask into the final result. We still recommend using the <span class="text-primary">
Python</span> language as he is quite handy.

```python
import json

import requests
import pandas as pd

baseUrl = 'http://121.40.210.46:8097'


def transform(data):
    resultList = []
    for mol in data['data']:
        if not mol['data']:
            # Invalid SMILES
            tmp = {'smiles': mol['smiles']}
        else:
            tmp = dict({'smiles': mol['smiles']})
            for _, admet in mol['data'].items():
                for endpoint in admet:
                    # endpoint is a dict
                    tmp[endpoint['name']] = endpoint['value']
        resultList.append(tmp)
    return pd.DataFrame(resultList).fillna('Invalid SMILES')


def divide_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


if __name__ == '__main__':
    api = '/api/admet'
    url = baseUrl + api
    param = {
        'SMILES': []
    }
    n = 1000
    smiles_list = ['CN1C2CCC1CC(OC(=O)c1cccn1C)C2', 'O=C(O)Nc1scnc1C(=O)Nc1nccs1',
                   'COc1ccc(C=C(F)C(=O)c2cc(OC)c(OC)c(OC)c2)cc1'] * 2500

    for _, sublist in enumerate(divide_list(smiles_list, n)):
        param['SMILES'] = sublist

        response = requests.post(url, json=param)

        if response.status_code == 200:  # If access is successful
            data = response.json()['data']
            # transform to csv file
            result = transform(data)
            result.to_csv('result' + str(_) + '.csv', index=False)
```