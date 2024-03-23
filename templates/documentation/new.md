# Evaluation Mode

---

## 1. Input Types

The <span class="blue-font">Evaluation Mode</span> allows the input of single molecule which aims
to <span class="blue-font">confirm the authenticity of biological campaign screening result</span>. Two
input types are provided:

1. By inputting single SMILES string;
2. By drawing molecule from editor below.

## 2. Single SMILES string

<span class="blue-font">Step1:</span> Access the Services page via Services-> Evaluation Mode in the navigation bar.

![](../static/docs/img/1.png ':class=content-img')

<span class="blue-font">Step2:</span> Select first entry for: Input a single strings and Enter the SMILES string in the
SMILES input box.

![](../static/docs/img/2.png ':class=content-img')

![](../static/docs/img/3.png ':class=content-img')

<span class="blue-font">Step3:</span> Select related evaluation methods for following application. Two types of
evaluation methods are provided: <span class="black-font">specific mechanism FH detection</span>
and <span class="black-font">collected credible FH screening rules application</span>. The prior type includes both
prediction models and substructure rules, which are constructed based on the special dataset collection.

![](../static/docs/img/4.png ':class=content-img')

![](../static/docs/img/5.png ':class=content-img')

<span class="blue-font">Mechanism:</span>

<ul class="page-content-ul">
    <li>
        1.<span class="black-font">Colloidal aggregators:</span> Compound aggregation tends to start when
        the concentration is above the CAC
        and end as aggregators form with radius of approximately 30-600 nm. The resulting colloidal
        aggregators would non-specifically bind to the surface of proteins, thus inducing local protein
        unfolding, which usually results in destabilization or denaturation of enzymes
    </li>
    <li>
        2.<span class="black-font">FLuc inhibitors:</span> Due to its unique catalysis mechanism, FLuc is
        widely used in a variety of HTS bioluminescence assays, especially in the assay which aims to study
        gene expression at the transcriptional level. However, the inhibition of Fluc by unexpected FLuc
        inhibitors would produce interference to HTS assays.
    </li>
    <li>
        3.<span class="black-font">Blue/Green fluorescence:</span> Fluorescence is the process by which a
        molecule, called fluorophore or fluorescent dye, absorbs a photon of light, exciting an electron to
        a higher energy state. Fluorophores have many applications, including as enzyme substrates, labels
        for biomolecules, cellular stains and environmental indicators. However, the appearance of
        fluorescent compound would produce interference to related HTS assays.
    </li>
    <li>
        4.<span class="black-font">Reactive compounds:</span> Chemical reactive compounds typically result
        in the chemical modification of reactive protein residues or, less frequently, the modification of
        nucleophilic assay reagents.
    </li>
    <li>
        5.<span class="black-font">Promiscuous compounds:</span> Promiscuous compounds refer to compounds
        that specifically bind to different macromolecular targets. These multiple interactions may include
        unintended targets, thus triggering adverse reactions and other safety issues.
    </li>
    <li>
        6.<span class="black-font">Other assay interferences:</span> Alpha-screen, FRET, TR-FRET, absorbance
        artifacts are included.
    </li>
</ul>

<span class="blue-font">Other FH Screening Rules:</span>

<ul class="page-content-ul">
    <li>
        1.<span class="black-font">PAINS:</span> frequent hitters, Alpha-screen artifacts and reactive
        compound; 480 substructures (J Med Chem 2010;53:2719–40)
    </li>
    <li>
        2.<span class="black-font">BMS:</span> undesirable, reactive compounds; 176 substructures (J Chem
        Inf Model 2006;46:1060–8)
    </li>
    <li>
        3.<span class="black-font">GST/GSH FH filter:</span> GST/GSH FHs; 34 substructures ( J Biomol Screen
        2016;21:596–607)
    </li>
    <li>
        4.<span class="black-font">His-tagged protein FH filter:</span> Ni2+ chelators; 19 substructures (J
        Biomol Screen 2014;19:715–26)
    </li>
    <li>
        5.<span class="black-font">ALARM NMR:</span> thiol reactive compounds; 75 substructures (J Am Chem
        Soc 2005;127:217–24)
    </li>
    <li>
        6.<span class="black-font">Luciferase inhibitor Rule:</span> FLuc inhibitors; 3 substructures (J
        Chem Inf Model 2018;58:933–42)
    </li>
    <li>
        7.<span class="black-font">Chelator Rule:</span> chelators; 55 substructures (ChemMedChem
        2010;5:195–9)
    </li>
    <li>
        8.<span class="black-font">NTD:</span> unwanted groups, reactive groups and possible HTS
        interferences; 105 substructures (ChemMedChem 2008;3:435–44)
    </li>
    <li>
        9.<span class="black-font">Potential electrophilic Rule:</span> reactive compounds; 119
        substructures (J Chem Inf Model 2012;52:2310–6)
    </li>
</ul>

<span class="blue-font">Step4:</span> Submit and get results.

![](../static/docs/img/6.png ':class=content-img')

## 3. Output Types

After submission, the backend will input the queried molecule into the selected models and (or) substructure rules. This
page will show a brief overview of the results obtained, including:

<span class="blue-font">1. Visualization block:</span>

<ul class="page-content-ul page-content-space">
    <li>
        1) The status of the queried molecule;
    </li>
    <li>
        2) The 2D, 3D structure graph of the queried molecule;
    </li>
    <li>
        3) The radar chart is provided for the better comprehension of the different mechanism prediction
        result, where Blue Label indicates that this type of mechanism is chosen for prediction, while Black
        Label indicates that this this type of mechanism is not chosen for prediction.
    </li>
</ul>

![](../static/docs/img/7.png ':class=content-img')

<span class="blue-font">2.Model Predictions and Substructure Screening blocks:</span>

<ul class="page-content-ul page-content-space">
    <li>
        1) <span class="black-font">Model Predictions block</span> provide the information about the
        prediction mechanism and corresponding
        GCN, GAT, AttentiveFP and average prediction scores;
    </li>
    <li>
        2) <span class="black-font">Substructure Screening block</span> provide the information about the
        prediction mechanism and corresponding flagged substructure name and graph.
    </li>
</ul>

![](../static/docs/img/8.png ':class=content-img')

<span class="blue-font">3.File Download:</span>

The csv (pdf) file shows the score the accepted status and fragment information for each mechanism for
each molecule (in this case, only one molecule), including SMILES, status, GCN score, GAT score,
AttentiveFP score, average score, and SMARTS for each rule.

![](../static/docs/img/9.png ':class=content-img')

## 4. Molecular Editor

<span class="blue-font">Step1:</span> Access the Services page via Services-> Evaluation Mode in the navigation bar.

![](../static/docs/img/1.png ':class=content-img')

<span class="blue-font">Step2:</span> Select second entry for: Drawing molecules through the molecule editor and draw a
molecule.

![](../static/docs/img/10.png ':class=content-img')

<span class="blue-font">Step3:</span> Select related evaluation methods for following application. Two
types of evaluation methods are provided: <span class="black-font">specific mechanism FH
detection</span> and <span class="black-font">collected credible FH screening rules application</span>. The prior type
includes both prediction models and substructure rules, which are constructed based on the special dataset collection.

![](../static/docs/img/11.png ':class=content-img')

![](../static/docs/img/12.png ':class=content-img')

<span class="blue-font">Mechanism:</span>

<ul class="page-content-ul">
    <li>
        1.<span class="black-font">Colloidal aggregators:</span> Compound aggregation tends to start when
        the concentration is above the CAC
        and end as aggregators form with radius of approximately 30-600 nm. The resulting colloidal
        aggregators would non-specifically bind to the surface of proteins, thus inducing local protein
        unfolding, which usually results in destabilization or denaturation of enzymes
    </li>
    <li>
        2.<span class="black-font">FLuc inhibitors:</span> Due to its unique catalysis mechanism, FLuc is
        widely used in a variety of HTS bioluminescence assays, especially in the assay which aims to study
        gene expression at the transcriptional level. However, the inhibition of Fluc by unexpected FLuc
        inhibitors would produce interference to HTS assays.
    </li>
    <li>
        3.<span class="black-font">Blue/Green fluorescence:</span> Fluorescence is the process by which a
        molecule, called fluorophore or fluorescent dye, absorbs a photon of light, exciting an electron to
        a higher energy state. Fluorophores have many applications, including as enzyme substrates, labels
        for biomolecules, cellular stains and environmental indicators. However, the appearance of
        fluorescent compound would produce interference to related HTS assays.
    </li>
    <li>
        4.<span class="black-font">Reactive compounds:</span> Chemical reactive compounds typically result
        in the chemical modification of reactive protein residues or, less frequently, the modification of
        nucleophilic assay reagents.
    </li>
    <li>
        5.<span class="black-font">Promiscuous compounds:</span> Promiscuous compounds refer to compounds
        that specifically bind to different macromolecular targets. These multiple interactions may include
        unintended targets, thus triggering adverse reactions and other safety issues.
    </li>
    <li>
        6.<span class="black-font">Other assay interferences:</span> Alpha-screen, FRET, TR-FRET, absorbance
        artifacts are included.
    </li>
</ul>

<span class="blue-font">Other FH Screening Rules:</span>

<ul class="page-content-ul">
    <li>
        1.<span class="black-font">PAINS:</span> frequent hitters, Alpha-screen artifacts and reactive
        compound; 480 substructures (J Med Chem 2010;53:2719–40)
    </li>
    <li>
        2.<span class="black-font">BMS:</span> undesirable, reactive compounds; 176 substructures (J Chem
        Inf Model 2006;46:1060–8)
    </li>
    <li>
        3.<span class="black-font">GST/GSH FH filter:</span> GST/GSH FHs; 34 substructures ( J Biomol Screen
        2016;21:596–607)
    </li>
    <li>
        4.<span class="black-font">His-tagged protein FH filter:</span> Ni2+ chelators; 19 substructures (J
        Biomol Screen 2014;19:715–26)
    </li>
    <li>
        5.<span class="black-font">ALARM NMR:</span> thiol reactive compounds; 75 substructures (J Am Chem
        Soc 2005;127:217–24)
    </li>
    <li>
        6.<span class="black-font">Luciferase inhibitor Rule:</span> FLuc inhibitors; 3 substructures (J
        Chem Inf Model 2018;58:933–42)
    </li>
    <li>
        7.<span class="black-font">Chelator Rule:</span> chelators; 55 substructures (ChemMedChem
        2010;5:195–9)
    </li>
    <li>
        8.<span class="black-font">NTD:</span> unwanted groups, reactive groups and possible HTS
        interferences; 105 substructures (ChemMedChem 2008;3:435–44)
    </li>
    <li>
        9.<span class="black-font">Potential electrophilic Rule:</span> reactive compounds; 119
        substructures (J Chem Inf Model 2012;52:2310–6)
    </li>
</ul>

<span class="blue-font">Step4:</span> Submit and get results.

![](../static/docs/img/13.png ':class=content-img')

## 5. Output Types

After submission, the backend will input the queried molecule into the selected models and (or)
substructure rules. This page will show a brief overview of the results obtained, including:

<span class="blue-font">1. Visualization block:</span>

<ul class="page-content-ul page-content-space">
    <li>
        1) The status of the queried molecule;
    </li>
    <li>
        2) The 2D, 3D structure graph of the queried molecule;
    </li>
    <li>
        3) The radar chart is provided for the better comprehension of the different mechanism prediction
        result, where Blue Label indicates that this type of mechanism is chosen for prediction, while Black
        Label indicates that this this type of mechanism is not chosen for prediction.
    </li>
</ul>

![](../static/docs/img/14.png ':class=content-img')

<span class="blue-font">2.Model Predictions and Substructure Screening blocks:</span>

<ul class="page-content-ul page-content-space">
    <li>
        1) <span class="black-font">Model Predictions block</span> provide the information about the
        prediction mechanism and corresponding
        GCN, GAT, AttentiveFP and average prediction scores;
    </li>
    <li>
        2) <span class="black-font">Substructure Screening block</span> provide the information about the
        prediction mechanism and corresponding flagged substructure name and graph.
    </li>
</ul>

![](../static/docs/img/15.png ':class=content-img')

<span class="blue-font">3.File Download:</span>

The csv (pdf) file shows the score the accepted status and fragment information for each mechanism for
each molecule (in this case, only one molecule), including SMILES, status, GCN score, GAT score,
AttentiveFP score, average score, and SMARTS for each rule.

![](../static/docs/img/16.png ':class=content-img')
