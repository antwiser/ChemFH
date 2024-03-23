# <svg t="1702051796815" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="12557" width="40" height="40"><path d="M781.0048 897.536H243.5584c-62.9248 0-113.9712-50.9952-113.9712-113.9712v-153.1904H894.976v153.1904c0 62.976-50.9952 113.9712-113.9712 113.9712z" fill="#ffa115" p-id="12558"></path><path d="M792.576 106.8544H231.9872c-73.4208 0-133.12 59.6992-133.12 133.12v560.5888c0 73.4208 59.6992 133.12 133.12 133.12H792.576c73.4208 0 133.12-59.6992 133.12-133.12V239.9744c0-73.3696-59.6992-133.12-133.12-133.12z m71.68 693.7088c0 39.5264-32.1536 71.68-71.68 71.68H231.9872c-39.5264 0-71.68-32.1536-71.68-71.68V239.9744c0-39.5264 32.1536-71.68 71.68-71.68H792.576c39.5264 0 71.68 32.1536 71.68 71.68v560.5888z" fill="#474A54" p-id="12559"></path><path d="M437.0432 306.6368H238.848c-16.9472 0-30.72 13.7728-30.72 30.72s13.7728 30.72 30.72 30.72h198.1952c16.9472 0 30.72-13.7728 30.72-30.72s-13.7728-30.72-30.72-30.72zM778.752 540.16h-198.144c-16.9472 0-30.72 13.7728-30.72 30.72s13.7728 30.72 30.72 30.72h198.144c16.9472 0 30.72-13.7728 30.72-30.72s-13.7728-30.72-30.72-30.72zM778.752 333.312h-68.352V264.96c0-16.9472-13.7728-30.72-30.72-30.72s-30.72 13.7728-30.72 30.72V333.312H580.608c-16.9472 0-30.72 13.7728-30.72 30.72s13.7728 30.72 30.72 30.72h68.352v68.352c0 16.9472 13.7728 30.72 30.72 30.72s30.72-13.7728 30.72-30.72V394.752h68.352c16.9472 0 30.72-13.7728 30.72-30.72s-13.7728-30.72-30.72-30.72zM434.4832 558.2848c-11.9808-11.9808-31.4368-11.9808-43.4688 0l-58.3168 58.3168-58.2656-58.3168c-11.9808-11.9808-31.4368-11.9808-43.4688 0-11.9808 11.9808-11.9808 31.4368 0 43.4688l58.3168 58.3168-58.2656 58.2656c-11.9808 11.9808-11.9808 31.4368 0 43.4688a30.67392 30.67392 0 0 0 43.4176 0l58.3168-58.3168 58.3168 58.3168a30.67392 30.67392 0 0 0 43.4176 0c11.9808-11.9808 11.9808-31.4368 0-43.4688l-58.3168-58.3168 58.3168-58.3168a30.62272 30.62272 0 0 0 0-43.4176zM778.752 697.9072h-198.144c-16.9472 0-30.72 13.7728-30.72 30.72s13.7728 30.72 30.72 30.72h198.144c16.9472 0 30.72-13.7728 30.72-30.72s-13.7728-30.72-30.72-30.72z" fill="#474A54" p-id="12560"></path></svg> ADMET Calculation

---

## 1. Calculate the ADMET properties of molecule(s).

> <span class="badge text-bg-primary">POST</span> <span color="grey">`http://121.40.210.46:8097`</span>/api/admet

<span style="color: #2878b5">Description:</span>

Calculate the ADMET properties of molecule(s), including molecular basic properties, physical chemistry, medicinal
chemistry, absorption, distribution, metabolism, excretion, and toxicity, among others.

> Each property includes computed values, detailed explanations, and provides decision results. In this context, 0
> represents <span class="text-success">'pass'</span>, -1 represents <span class="text-warning">'warning'</span>, and 1
> represents <span class="text-danger">'danger'</span>. Any other numbers indicate that the property has no decision.
> Additionally, for properties related to toxicity groups or substructure rule matches, the decision value may be an
> array
> of SVG strings, representing highlighted images of molecules matching the endpoint. These 'values' are displayed as
> two-dimensional arrays, where each element represents the indices of highlighted matching atoms in the molecule.

<span style="color: #2878b5">Query Parameters & Example:</span>

| Key     | Value                                                 | Type               | Description                                                                    |
|---------|-------------------------------------------------------|--------------------|--------------------------------------------------------------------------------|
| SMILES  | ["molecule", "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"] | string or [string] | Molecular SMILES string, which is commonly provided by most chemical toolkits. |
| feature | false                                                 | boolean            | Whether to use models with descriptors                                         |

<span style="color: #2878b5">Result & Example:</span>

Return all computable properties of the molecule, 'structure' represents the SVG string of the molecular structure,
and 'taskid' serves as the unique identifier of the task, which can be used to retrieve the result file later.

Content-Type: `application/json,application/xml`

```json
{
    "status": "success",
    "code": 200,
    "data": {
        // The unique identifier for a submitted task.
        "taskid": "tmpvvsp7s291702724373", 
        // Number of all molecules.
        "number of all molecules": 1,
        // Number of valid molecules.
        "number of valid molecules": 0,
        // Response data.
        "data": [
            {
                "smiles": "molecule",
                "MW": "Invalid Molecule",
                "Vol": "Invalid Molecule",
                "Dense": "Invalid Molecule",
                "nHA": "Invalid Molecule",
                "nHD": "Invalid Molecule",
                "TPSA": "Invalid Molecule"
                // ...
            },
            {
                "smiles": "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2",
                "MW": 318.1,
                "Vol": 316.59731852400796,
                "Dense": 1.0047463493468538,
                "nHA": 5,
                "nHD": 0,
                "TPSA": 80.05
                // ...
            }
        ],
        // Explanation of each key.
        "explanation": {
            "Aggregators": "Category 0: non-colloidal aggregators; Category 1: colloidal aggregators. The output value is the probability of being colloidal aggregators, within the range of 0 to 1.",
            "Fluc": "Category 0: non-fLuc inhibitors; Category 1: fLuc inhibitors. The output value is the probability of being fLuc inhibitors, within the range of 0 to 1.",
            "Blue_fluorescence": "Category 0: non-blue fluorescence; Category 1: blue fluorescence. The output value is the probability of being blue fluorescence, within the range of 0 to 1.",
            "Green_fluorescence": "Category 0: non-green fluorescence; Category 1: green fluorescence. The output value is the probability of being green fluorescence, within the range of 0 to 1.",
            "Reactive": "Category 0: non-reactive compound; Category 1: reactive compound. The output value is the probability of being reactive compound, within the range of 0 to 1.",
            "Other_assay_interference": "",
            // ...
        }
    }
}
```

|     |                       |                                                               |
|-----|-----------------------|---------------------------------------------------------------|
| 200 | OK                    | The service call has completed successfully.                  |
| 500 | Internal Server Error | An unexpected error has happened while resolving the request. |

## 2. Retrieve a CSV file of the calculation results for molecule(s).

> <span class="badge text-bg-primary">POST</span> <span color="grey">`http://121.40.210.46:8097`</span>/api/admetCSV

<span style="color: #2878b5">Description:</span>

Retrieve the result file of the calculation by using the task ID. Each line in the result file represents a molecule,
and if the SMILES of a molecule is invalid, the result will display "Invalid SMILES".

<span style="color: #2878b5">Query Parameters & Example:</span>

| Key    | Value                 | Type   | Description                                             |
|--------|-----------------------|--------|---------------------------------------------------------|
| taskId | tmpva_lbtw11700216887 | string | The 'taskId' in the returned result of the request 2.1. |

<span style="color: #2878b5">Result & Example:</span>

CSV file of ADMET results.

Content-Type: `text/csv`

<div class="text-center">
<a href="../static/apis/files/results.csv" download="result.csv">
    <button type="button" class="btn btn-success"><i class="bi bi-cloud-arrow-down-fill"></i> Download Example</button>
</a></div>

