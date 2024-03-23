# <svg t="1702051514762" class="icon" viewBox="0 0 1032 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2342" width="40" height="40"><path d="M514.285714 1024a514.285714 514.285714 0 0 1-514.285714-509.714286 514.285714 514.285714 0 0 1 1028.571429 0 514.285714 514.285714 0 0 1-514.285715 509.714286zM514.285714 61.142857a457.142857 457.142857 0 0 0-457.142857 453.142857 457.142857 457.142857 0 0 0 914.285714 0A457.142857 457.142857 0 0 0 514.285714 61.142857z" fill="#0F3567" p-id="2343"></path><path d="M510.285714 486.857143a28.571429 28.571429 0 0 1-18.857143-6.857143L356.571429 359.428571a28.571429 28.571429 0 0 1-8-29.714285l36-120.571429a28.571429 28.571429 0 0 1 27.428571-20h189.714286a28 28 0 0 1 26.857143 20.571429l36 124.571428a28 28 0 0 1-8 28.571429l-127.428572 114.285714a28.571429 28.571429 0 0 1-18.857143 9.714286zM408.571429 328.571429l101.714285 91.428571 94.857143-86.857143-25.142857-86.857143H433.142857z" fill="#0F3567" p-id="2344"></path><path d="M506.857143 909.142857a28.571429 28.571429 0 0 1-21.714286-9.714286L298.285714 684a28 28 0 0 1-4.571428-30.285714l91.428571-206.857143 52 22.857143-84 190.285714 153.714286 177.142857 148.571428-171.428571-65.714285-198.285715 54.285714-18.285714 70.857143 211.428572a28 28 0 0 1-5.714286 27.428571l-180.571428 211.428571a28.571429 28.571429 0 0 1-21.714286 9.714286zM675.142857 442.4l53.257143-20.742857 41.885714 107.542857-53.257143 20.742857z" fill="#0F3567" p-id="2345"></path></svg> Introduction

---

?> This document aims to elucidate how ChemFH's API service is structured, with a variety of usage cases as
illustrations, to help new users learn how the service works and how to construct the URLs that are the interface to
this service.

> **USAGE POLICY:** Please note that ChemFH API is not designed for processing very large volumes (millions) of
> requests at the same time. We suggest that any script or application refrain from exceeding 5 requests per second to
> prevent overloading the ChemFH servers. If you have a large data set that you need to compute with, please contact
> us for assistance in optimizing your task, as there are likely more efficient ways to approach such bulk queries.

![ChemFH Framework](../static/apis/img/api_overview.png ':class=content-img')

<center><strong>Figure 1.</strong> Overview scheme of ChemFH API</center>

### How ChemFH Works

---


The conceptual framework of ChemFH is primarily composed of three-part request: **1) input**—cleaning individual
molecules or batches of molecules; **2) operation**—calculating the ADMET properties of these molecules using deep
learning models; and **3) output**—determining which result files to return. The beauty of this design is that all the
underlying functionalities are modular, allowing for flexible combinations.

![](../static/apis/img/dev-overview.png)

For example, this service supports input of chemical structure by SMILES. It supports output of chemical structure as
images in SVG format. You can combine these two into a visualization request for a SMILES string. - In this case, you
can obtain the SVG image of the molecule by sending a <span class="badge text-bg-primary">POST</span> request with the
SMILES string to ```http://121.40.210.46:8097/api/molsvg``` using the following request format:

```json
{
    "SMILES": "CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N"
}
```

We will get the output:

![](../static/apis/img/example.svg)

The possibilities are nearly endless, and more importantly, the action of the service is simple to understand from the
URL alone, without needing any extra programming, parsing, etc. And at the same time, more complex data handling is
available for programmers who want rigorous schema-based XML communications, or who want to use JSON data to embed
functionality in a web page via JavaScript.
