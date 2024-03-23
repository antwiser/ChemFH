# How ADMETlab 3.0 Works

---

The conceptual framework of ADMETlab 3.0 is primarily composed of three-part request: **1) input**—cleaning individual molecules or batches of molecules; **2) operation**—calculating the ADMET properties of these molecules using deep learning models; and **3) output**—determining which result files to return. The beauty of this design is that all the underlying functionalities are modular, allowing for flexible combinations.

![](../static/apis/img/dev-overview.png)

For example, this service supports input of chemical structure by SMILES. It supports output of chemical structure as images in SVG format. You can combine these two into a visualization request for a SMILES string. - In this case, you can obtain the SVG image of the molecule by sending a <span class="badge text-bg-primary">POST</span> request with the SMILES string to ```http://121.40.210.46:8098/api/molsvg``` using the following request format:

```json
{ "SMILES": "CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N" }
```

We will get the output:

![](../static/apis/img/example.svg)

The possibilities are nearly endless, and more importantly, the action of the service is simple to understand from the URL alone, without needing any extra programming, parsing, etc. And at the same time, more complex data handling is available for programmers who want rigorous schema-based XML communications, or who want to use JSON data to embed functionality in a web page via JavaScript.

