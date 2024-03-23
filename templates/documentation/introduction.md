# Introduction

---

## Evaluation Mode

The <span class="blue-font">Evaluation Mode</span> allows the input of <span class="black-font">single molecule</span>
which aims to confirm the authenticity of biological campaign screening result. Based on the combination of credible
prediction models and useful substructure rules, it is believed that through Evaluation Mode, it is able to
comprehensively evaluate the risk and potential of the queried molecule.

## Screening Mode

The <span class="blue-font">Screening Mode</span> is able to assist hit identification by
screening <span class="black-font">molecular dataset</span> prior to detect potential frequent hitters.
Such function makes it possible to pre-screen molecules with undesirable interfering features, thus increasing the
efficiency of drug design, improving the credibility of experimental result and decreasing unnecessary cost.

## Performance of Prediction Models

<div class="table-wrapper">
<table>
    <tr class="blue-header">
        <td rowspan="2" class="left aligned middle aligned">
            <p>
                Category
            </p>
        </td>
        <td colspan="3" class="center aligned">
            <p>
                Data
            </p>
        </td>
        <td colspan="3" class="center aligned">
            <p>
                Model (AUC)
            </p>
        </td>
    </tr>
    <tr class="blue-header">
        <td class="center aligned">
            <p>
                Positive
            </p>
        </td>
        <td class="center aligned">
            <p>
                Negative
            </p>
        </td>
        <td class="center aligned">
            <p>
                Total
            </p>
        </td>
        <td class="center aligned">
            <p>
                GCN
            </p>
        </td>
        <td class="center aligned">
            <p>
                GAT
            </p>
        </td>
        <td class="center aligned">
            <p>
                AttentiveFP
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                Agg
            </p>
        </td>
        <td>
            <p>
                16783
            </p>
        </td>
        <td>
            <p>
                35114
            </p>
        </td>
        <td>
            <p>
                51897
            </p>
        </td>
        <td>
            <p>
                0.937
            </p>
        </td>
        <td>
            <p>
                0.937
            </p>
        </td>
        <td>
            <p>
                0.937
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                Blue FLuo
            </p>
        </td>
        <td>
            <p>
                4871
            </p>
        </td>
        <td>
            <p>
                38659
            </p>
        </td>
        <td>
            <p>
                43530
            </p>
        </td>
        <td>
            <p>
                0.936
            </p>
        </td>
        <td>
            <p>
                0.938
            </p>
        </td>
        <td>
            <p>
                0.926
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                GreenFLuo
            </p>
        </td>
        <td>
            <p>
                8319
            </p>
        </td>
        <td>
            <p>
                31352
            </p>
        </td>
        <td>
            <p>
                39671
            </p>
        </td>
        <td>
            <p>
                0.728
            </p>
        </td>
        <td>
            <p>
                0.723
            </p>
        </td>
        <td>
            <p>
                0.720
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                Other assay interference
            </p>
        </td>
        <td>
            <p>
                1544
            </p>
        </td>
        <td>
            <p>
                2159
            </p>
        </td>
        <td>
            <p>
                3703
            </p>
        </td>
        <td>
            <p>
                0.779
            </p>
        </td>
        <td>
            <p>
                0.777
            </p>
        </td>
        <td>
            <p>
                0.743
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                FLuc inhibitor
            </p>
        </td>
        <td>
            <p>
                12703
            </p>
        </td>
        <td>
            <p>
                121345
            </p>
        </td>
        <td>
            <p>
                133418
            </p>
        </td>
        <td>
            <p>
                0.969
            </p>
        </td>
        <td>
            <p>
                0.967
            </p>
        </td>
        <td>
            <p>
                0.910
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                Reactive compound
            </p>
        </td>
        <td>
            <p>
                820
            </p>
        </td>
        <td>
            <p>
                162069
            </p>
        </td>
        <td>
            <p>
                162889
            </p>
        </td>
        <td>
            <p>
                0.992
            </p>
        </td>
        <td>
            <p>
                0.993
            </p>
        </td>
        <td>
            <p>
                0.909
            </p>
        </td>
    </tr>
    <tr>
        <td class="left aligned">
            <p>
                Promiscuous compound
            </p>
        </td>
        <td>
            <p>
                7518
            </p>
        </td>
        <td>
            <p>
                334171
            </p>
        </td>
        <td>
            <p>
                341689
            </p>
        </td>
        <td>
            <p>
                0.921
            </p>
        </td>
        <td>
            <p>
                0.930
            </p>
        </td>
        <td>
            <p>
                0.905
            </p>
        </td>
    </tr>
</table>
</div>

## Performance of Mechanism Substructure Rules

<div class="table-wrapper">
<table>
<tr class="blue-header">
    <td class="left aligned middle aligned">
        <p>
            Category
        </p>
    </td>
    <td>
        <p>
            Number of substructures
        </p>
    </td>
    <td>
        <p>
            Flagged molecules
        </p>
    </td>
    <td>
        <p>
            Flagged FHs
        </p>
    </td>
    <td>
        <p>
            Accuracy
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            Agg
        </p>
    </td>
    <td>
        <p>
            26
        </p>
    </td>
    <td>
        <p>
            2651
        </p>
    </td>
    <td>
        <p>
            2265
        </p>
    </td>
    <td>
        <p>
            0.854
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            Blue FLuo
        </p>
    </td>
    <td>
        <p>
            26
        </p>
    </td>
    <td>
        <p>
            1208
        </p>
    </td>
    <td>
        <p>
            1172
        </p>
    </td>
    <td>
        <p>
            0.970
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            GreenFLuo
        </p>
    </td>
    <td>
        <p>
            16
        </p>
    </td>
    <td>
        <p>
            274
        </p>
    </td>
    <td>
        <p>
            260
        </p>
    </td>
    <td>
        <p>
            0.949
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            Other assay interference
        </p>
    </td>
    <td>
        <p>
            7
        </p>
    </td>
    <td>
        <p>
            209
        </p>
    </td>
    <td>
        <p>
            200
        </p>
    </td>
    <td>
        <p>
            0.957
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            FLuc inhibitor
        </p>
    </td>
    <td>
        <p>
            18
        </p>
    </td>
    <td>
        <p>
            1638
        </p>
    </td>
    <td>
        <p>
            1448
        </p>
    </td>
    <td>
        <p>
            0.884
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            Reactive compound
        </p>
    </td>
    <td>
        <p>
            3
        </p>
    </td>
    <td>
        <p>
            30
        </p>
    </td>
    <td>
        <p>
            23
        </p>
    </td>
    <td>
        <p>
            0.767
        </p>
    </td>
</tr>
<tr>
    <td class="left aligned">
        <p>
            Promiscuous compound
        </p>
    </td>
    <td>
        <p>
            6
        </p>
    </td>
    <td>
        <p>
            733
        </p>
    </td>
    <td>
        <p>
            704
        </p>
    </td>
    <td>
        <p>
            0.960
        </p>
    </td>
</tr>
</table>
</div>

## Information about Other FH Screening Rules

<div class="table-wrapper">
<table>
<tr class="blue-header">
    <td width="20%">
        <p>
            Name
        </p>
    </td>
    <td width="25%">
        <p>
            Number of substructures
        </p>
    </td>
    <td width="30%">
        <p>
            Screening Scope
        </p>
    </td>
    <td>
        <p>
            Source
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            PAINS
        </p>
    </td>
    <td>
        <p class="text-center">
            480
        </p>
    </td>
    <td>
        <p>
            Alpha-screen artifacts and reactive compound
        </p>
    </td>
    <td>
        <p>
            J Med Chem 2010;53:2719–40
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            BMS
        </p>
    </td>
    <td>
        <p class="text-center">
            176
        </p>
    </td>
    <td>
        <p>
            Undesirable, reactive compounds
        </p>
    </td>
    <td>
        <p>
            J Chem Inf Model 2006;46:1060–8
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            GST/GSH FH filter
        </p>
    </td>
    <td>
        <p class="text-center">
            34
        </p>
    </td>
    <td>
        <p>
            GST/GSH FHs
        </p>
    </td>
    <td>
        <p>
            J Biomol Screen 2016;21:596–607
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            His-tagged protein FH filter
        </p>
    </td>
    <td>
        <p class="text-center">
            19
        </p>
    </td>
    <td>
        <p>
            Ni2+&nbsp;chelators
        </p>
    </td>
    <td>
        <p>
            J Biomol Screen 2014;19:715–26
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            ALARM NMR
        </p>
    </td>
    <td>
        <p class="text-center">
            75
        </p>
    </td>
    <td>
        <p>
            Thiol reactive compounds
        </p>
    </td>
    <td>
        <p>
            J Am Chem Soc 2005;127:217–24
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            Luciferase inhibitor Rule
        </p>
    </td>
    <td>
        <p class="text-center">
            3
        </p>
    </td>
    <td>
        <p>
            FLuc inhibitors
        </p>
    </td>
    <td>
        <p>
            J Chem Inf Model 2018;58:933–42
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            Chelator Rule
        </p>
    </td>
    <td>
        <p class="text-center">
            55
        </p>
    </td>
    <td>
        <p>
            Chelators
        </p>
    </td>
    <td>
        <p>
            ChemMedChem 2010;5:195–9
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            NTD
        </p>
    </td>
    <td>
        <p class="text-center">
            105
        </p>
    </td>
    <td>
        <p>
            Unwanted groups, reactive groups and possible HTS interferences
        </p>
    </td>
    <td>
        <p>
            ChemMedChem 2008;3:435–44
        </p>
    </td>
</tr>
<tr>
    <td>
        <p>
            Potential electrophilic Rule
        </p>
    </td>
    <td>
        <p class="text-center">
            119
        </p>
    </td>
    <td>
        <p>
            Reactive compounds
        </p>
    </td>
    <td>
        <p>
            J Chem Inf Model 2012;52:2310–6
        </p>
    </td>
</tr>
</table>
</div>

## Development Environment

<div class="table-wrapper">
<table>
    <tr class="blue-header">
        <td>
            <p>
                Library
            </p>
        </td>
        <td>
            <p>
                Version
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                RDkit
            </p>
        </td>
        <td>
            <p>
                2019.03.1
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                Django
            </p>
        </td>
        <td>
            <p>
                2.2
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                DGL
            </p>
        </td>
        <td>
            <p>
                0.5.2
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                DGL-LifeSci
            </p>
        </td>
        <td>
            <p>
                0.2.5
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                Pytorch
            </p>
        </td>
        <td>
            <p>
                1.6.0
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                Torchvision
            </p>
        </td>
        <td>
            <p>
                0.7.0
            </p>
        </td>
    </tr>
    <tr>
        <td>
            <p>
                pyecharts
            </p>
        </td>
        <td>
            <p>
                1.8.1
            </p>
        </td>
    </tr>
</table>
</div>