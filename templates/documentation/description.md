# Program Description

---

## 1. Model Datasets

<center><strong>Table 1.</strong> Data information of seven frequent hitter mechanism models.</center>

<div class="table-wrapper">
    <table>
        <thead>
        <tr>
            <th width="20%">Dataset</th>
            <th>Total (Positive/Negative)</th>
            <th>Training set (Positive/Negative)</th>
            <th>Test set (Positive/Negative)</th>
            <th>Validation set (Positive/Negative)</th>
        </tr>
        </thead>
        <tbody>
        <tr>
            <td>Colloidal aggregators</td>
            <td>51952 (16820/35132)</td>
            <td>41562 (13456/28106)</td>
            <td>5195 (1682/3513)</td>
            <td>5195 (1682/3513)</td>
        </tr>
        <tr>
            <td>FLuc inhibitors</td>
            <td>133452 (12092/121360)</td>
            <td>106762 (9674/97088)</td>
            <td>13345 (1209/12136)</td>
            <td>13345 (1209/12136)</td>
        </tr>
        <tr>
            <td>Blue fluorescent</td>
            <td>43597 (4871/38726)</td>
            <td>34877 (3897/30980)</td>
            <td>4360 (487/3873)</td>
            <td>4360 (487/3873)</td>
        </tr>
        <tr>
            <td>Green fluorescent</td>
            <td>47288 (8441/38847)</td>
            <td>37830 (6753/31077)</td>
            <td>4729 (844/3885)</td>
            <td>4729 (844/3885)</td>
        </tr>
        <tr>
            <td>Reactive compounds</td>
            <td>169761 (6056/163705)</td>
            <td>135809 (4844/130965)</td>
            <td>16976 (606/16370)</td>
            <td>16976 (606/16370)</td>
        </tr>
        <tr>
            <td>Promiscuous compounds</td>
            <td>373598 (8089/365509)</td>
            <td>298878 (6471/292407)</td>
            <td>37360 (809/36551)</td>
            <td>37360 (809/36551)</td>
        </tr>
        <tr>
            <td>Assay interference</td>
            <td>3743 (1544/2199)</td>
            <td>2995 (1236/1759)</td>
            <td>374 (154/220)</td>
            <td>374 (154/220)</td>
        </tr>
        </tbody>
    </table>
</div>

## 2. DMPNN Framework

Yang et al.(1). introduced an open-source Python package called Chemprop designed for implementing DMPNN models. This
package offers a robust and efficient solution tailored for molecular property prediction tasks and has garnered
extensive utilization in fields such as drug discovery and materials science(2). Within the DMPNN framework,
there exist two distinct stages: message passing and readout. Here, a graph G serves as an illustrative example,
representing node (atom) features as $x_v$ and edge (bond) features as $e_{vw}$. Initially, the edge hidden state $h_{v
w}^0$ is initialized using eq 1. This equation concatenate atom and bond features by passing them through the learned
matrix
$W_I$ and applying the rectified linear unit (RELU) activation function. This initialization defines the edge hidden
state,
which undergoes subsequent updates during the message passing process. In eq 1, τ represent the RELU activation
function, $W_I$ denote a learned matrix, and $cat()$ signifies a basic concatenation operation.

$$
h_{v w}^0=\tau\left(W_i \operatorname{cat}\left(x_v, e_{v w}\right)\right) \textcolor{red}{(1)}
$$

The first phase of message passing computes interactions between atom $v$ and atom $w$. This is achieved by summing the
hidden states of all bonds connected to atom v while excluding the hidden state of bonds from atom $w$, as described by
$m_{vw}^{(t+1)}$ in eq 2. This step captures information about the neighboring atoms and their connections to
individual atoms
within vw. In eq 2, $h_{kv}^t∈R^D$ represents bond features at layer $t∈{1,2,…,T}, k∈{N(v)\\w}$ denotes the set of nodes
connected to $v$ excluding $w$.

$$
m_{v w}^{t+1}=\sum_{k \in\{N(v) \backslash w\}} h_{k v}^t \textcolor{red}{(2)}
$$

Following this, a new hidden message at depth 1 is created by summing the product of the initial hidden state and the
learned matrix $W_m$ with the message. This resultant output undergoes further processing using the activation function
$τ$, denoted as $h_{vw}^{(t+1)}$ in eq 3.

$$ h_{v w}^{t+1}=\tau\left(h_{v w}^0+W_m m_{v w}^{t+1}\right) \textcolor{red}{(3)} $$

In the final message passing layer (at $t=T$), the updated hidden states $h_{vw}^T$ are summed to create the ultimate
message
for each atom, as described in eq 4. This action aggregates information about all neighboring atoms and their
relationships into the final message for each atom.

$$ m_v=\sum_{w \in N(v)} h_{v w}^T \textcolor{red}{(4)} $$

The hidden state $h_v$ for each atom is derived by concatenating the initial atom features with the message vector,
as indicated in eq 5.

$$
h_v=\tau\left(W_i \operatorname{cat}\left(x_v, m_v\right)\right) \textcolor{red}{(5)}
$$

Finally, employing eq 6, the hidden states $h_v$ of each atom are summed to generate a molecular feature vector. This
step
aggregate information from all atoms in the molecule into a unified molecular feature vector, facilitating property
prediction. It encompasses both structural and attribute information of the entire molecule, offering a comprehensive
representation for further property prediction.

$$
h=\sum_{v \in N(v)} h_v \textcolor{red}{(6)}
$$

In the DMPNN-Des model, preceding the readout phase, vector $h$ is concatenated with descriptor vectors and collectively
processed using a fully connected feedforward neural network to predict properties. The algorithm is implemented using
the open-source Chemprop package(2). Regarding the raw datasets, we trained each dataset using Chemprop, employing
random segmentation ratios [0.8, 0.1, 0.1] for training, testing, and evaluation. The batch process is iterated five
times, and the average RMSE value and variance from these iterations are calculated to evaluate the robustness of the
model.

## 3. Model performance

<center><strong>Table 2.</strong> Performance of MPNN-Des, DMPNN-FP, DMPNN on the Test Set.</center>

<div class="table-wrapper">
<table>
    <thead>
    <tr>
        <th>Dataset</th>
        <th>Model</th>
        <th>BA</th>
        <th>ACC</th>
        <th>AUC</th>
        <th>SP</th>
        <th>SE</th>
        <th>MCC</th>
    </tr>
    </thead>
    <tbody>
    <tr>
        <td rowspan="3">Aggregators</td>
        <td>DMPNN-Des</td>
        <td>0.871±0.005</td>
        <td>0.891±0.004</td>
        <td>0.938±0.004</td>
        <td>0.928±0.015</td>
        <td>0.813±0.022</td>
        <td>0.749±0.006</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.862±0.004</td>
        <td>0.886±0.004</td>
        <td>0.934±0.003</td>
        <td>0.929±0.013</td>
        <td>0.796±0.019</td>
        <td>0.737±0.007</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.861±0.005</td>
        <td>0.882±0.004</td>
        <td>0.933±0.003</td>
        <td>0.921±0.013</td>
        <td>0.801±0.019</td>
        <td>0.729±0.009</td>
    </tr>
    <tr>
        <td rowspan="3">FLuc inhibitors</td>
        <td>DMPNN-Des</td>
        <td>0.908±0.006</td>
        <td>0.971±0.001</td>
        <td>0.979±0.002</td>
        <td>0.985±0.002</td>
        <td>0.831±0.013</td>
        <td>0.824±0.008</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.906±0.008</td>
        <td>0.971±0.002</td>
        <td>0.978±0.003</td>
        <td>0.986±0.001</td>
        <td>0.827±0.016</td>
        <td>0.823±0.011</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.901±0.009</td>
        <td>0.971±0.001</td>
        <td>0.979±0.002</td>
        <td>0.986±0.002</td>
        <td>0.817±0.019</td>
        <td>0.820±0.007</td>
    </tr>
    <tr>
        <td rowspan="3">Blue fluorescent compounds</td>
        <td>DMPNN-Des</td>
        <td>0.874±0.008</td>
        <td>0.953±0.003</td>
        <td>0.947±0.004</td>
        <td>0.976±0.005</td>
        <td>0.773±0.019</td>
        <td>0.761±0.015</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.865±0.012</td>
        <td>0.952±0.003</td>
        <td>0.941±0.007</td>
        <td>0.977±0.004</td>
        <td>0.753±0.026</td>
        <td>0.752±0.016</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.864±0.014</td>
        <td>0.953±0.004</td>
        <td>0.946±0.005</td>
        <td>0.979±0.005</td>
        <td>0.748±0.03</td>
        <td>0.757±0.019</td>
    </tr>
    <tr>
        <td rowspan="3">Green fluorescent compounds</td>
        <td>DMPNN-Des</td>
        <td>0.659±0.008</td>
        <td>0.751±0.020</td>
        <td>0.729±0.007</td>
        <td>0.801±0.033</td>
        <td>0.517±0.041</td>
        <td>0.281±0.014</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.654±0.006</td>
        <td>0.747±0.040</td>
        <td>0.722±0.006</td>
        <td>0.799±0.062</td>
        <td>0.509±0.064</td>
        <td>0.275±0.021</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.661±0.006</td>
        <td>0.744±0.042</td>
        <td>0.727±0.008</td>
        <td>0.790±0.064</td>
        <td>0.532±0.064</td>
        <td>0.283±0.026</td>
    </tr>
    <tr>
        <td rowspan="3">Reactive compounds</td>
        <td>DMPNN-Des</td>
        <td>0.972±0.006</td>
        <td>0.972±0.005</td>
        <td>0.992±0.003</td>
        <td>0.976±0.005</td>
        <td>0.967±0.010</td>
        <td>0.944±0.010</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.970±0.003</td>
        <td>0.971±0.003</td>
        <td>0.993±0.002</td>
        <td>0.977±0.006</td>
        <td>0.964±0.008</td>
        <td>0.941±0.006</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.973±0.004</td>
        <td>0.973±0.004</td>
        <td>0.991±0.002</td>
        <td>0.974±0.007</td>
        <td>0.972±0.008</td>
        <td>0.946±0.008</td>
    </tr>
    <tr>
        <td rowspan="3">Promiscuous compounds</td>
        <td>DMPNN-Des</td>
        <td>0.877±0.007</td>
        <td>0.878±0.008</td>
        <td>0.936±0.007</td>
        <td>0.892±0.019</td>
        <td>0.862±0.013</td>
        <td>0.754±0.016</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.874±0.008</td>
        <td>0.875±0.009</td>
        <td>0.935±0.006</td>
        <td>0.881±0.020</td>
        <td>0.867±0.014</td>
        <td>0.748±0.018</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.880±0.006</td>
        <td>0.882±0.006</td>
        <td>0.937±0.005</td>
        <td>0.890±0.015</td>
        <td>0.871±0.017</td>
        <td>0.761±0.012</td>
    </tr>
    <tr>
        <td rowspan="3">Assay interference</td>
        <td>DMPNN-Des</td>
        <td>0.777±0.032</td>
        <td>0.770±0.037</td>
        <td>0.844±0.027</td>
        <td>0.736±0.066</td>
        <td>0.818±0.036</td>
        <td>0.548±0.063</td>
    </tr>
    <tr>
        <td>DMPNN-FP</td>
        <td>0.764±0.015</td>
        <td>0.750±0.021</td>
        <td>0.832±0.014</td>
        <td>0.684±0.06</td>
        <td>0.845±0.048</td>
        <td>0.524±0.028</td>
    </tr>
    <tr>
        <td>DMPNN</td>
        <td>0.774±0.035</td>
        <td>0.764±0.046</td>
        <td>0.846±0.027</td>
        <td>0.717±0.102</td>
        <td>0.831±0.044</td>
        <td>0.545±0.065</td>
    </tr>
    </tbody>
</table>
</div>

## 4. Uncertainty Estimation

In the ChemFH model, uncertainty is estimated using the Monte Carlo dropout approach(3). This method considers
dropout in deep neural networks as an approximate Bayesian inference of deep Gaussian processes. Specifically, the
dropout technique involves applying dropout before each layer during training and maintaining dropout activation during
the inference process. This allows for the generation of prediction distributions using different random masks,
approximating the posterior of deep Gaussian processes. The variance of this distribution serves as an
estimation of predictive uncertainty(4-6).

During experiments with uncertainty quantification, an ensemble of models generates a prediction distribution for each
molecule using a dropout-enabled network with a sample size of 10. The probability of 0.1 to use for Monte Carlo dropout
uncertainty estimation. Let $y_t$ represent the prediction from a single model within the ensemble, which contains
$T = 10$
models. For a query sample $x_i$, the prediction $\hat{y}$ is represented as the means of all predictions and the
uncertainty of
this sample U(x) can be provided by the variance $σ_t^2$ of the prediction distribution.

$$ U(x)=\sigma^2=\frac{1}{T-1} \sum_{t=1}^T\left(y_t-\hat{y}\right)^2 $$

We employed the method proposed by Dolezal et al.(7). to determine the optimal uncertainty threshold value, $θ$. This
method establishes an uncertainty threshold, where predictions below this threshold are more likely to be correct than
those with higher levels of uncertainty. To find the uncertainty threshold that optimally separates predictions into
likely-correct (high-confidence) and likely-incorrect (low-confidence), we calculated the sensitivity and specificity
for misprediction for all possible uncertainty thresholds. The corresponding Youden’s index (J) for each uncertainty
threshold $θ_i$ is the calculated as

$$
J_i=S e_i+S p_i-1
$$

The optimal uncertainty threshold $θ$ is the defined as the threshold which maximized the Youden's index:

$$
\theta=\operatorname{argmax}_i
$$

The single threshold is then used for all predicted made by the model. We take a binary approach to confidence
using the uncertainty threshold, with confidence of the classification model defined as

$$
C(x)= \begin{cases}\text { high }- \text { confidence } & \sigma^2(x)<\theta \\ \text { low }- \text { confidence } &
\sigma^2(x) \geq \theta\end{cases}
$$

In other words, prediction uncertainty exceeding this value designates the model's prediction as
<span class='low-content'>"Low-confidence"</span>,
while
prediction uncertainty below this threshold indicates <span class='high-content'>"High-confidence"</span> in the model's
prediction. This
threshold will
be used to assess the reliability of prediction in classification tasks within the ChemFH models.

## 5. Running Time

<center><strong>Table 3.</strong> Runtime analysis in seconds for submissions of 1 to 10000 molecules for ChemFH.</center>

<div class="table-wrapper">
<table>
    <thead>
    <tr>
        <th>Number of molecules</th>
        <th>ChemFH</th>
    </tr>
    </thead>
    <tbody>
    <tr>
        <td>1</td>
        <td>2.03</td>
    </tr>
    <tr>
        <td>10</td>
        <td>3.23</td>
    </tr>
    <tr>
        <td>100</td>
        <td>4.12</td>
    </tr>
    <tr>
        <td>500</td>
        <td>6.44</td>
    </tr>
    <tr>
        <td>1,000</td>
        <td>8.51</td>
    </tr>
    <tr>
        <td>5,000</td>
        <td>26.52</td>
    </tr>
    <tr>
        <td>10,000</td>
        <td>50.69</td>
    </tr>
    </tbody>
</table>
</div>

## 6. Browser Compatibility

<center><strong>Table 4.</strong> The information about the development environment of ChemFH.
</center>

<div class="table-wrapper">
    <table>
        <thead>
        <tr>
            <th>Name</th>
            <th>Version</th>
            <th class="text-center">Source</th>
        </tr>
        </thead>
        <tbody>
        <tr>
            <td>RDKit</td>
            <td>2023.9.4</td>
            <td class="text-center"><a href="https://www.rdkit.org/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>django</td>
            <td>4.1</td>
            <td class="text-center"><a href="https://www.djangoproject.com/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>django-ninja</td>
            <td>1.1.0</td>
            <td class="text-center"><a href="https://django-ninja.dev/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>dgl</td>
            <td>1.1.2</td>
            <td class="text-center"><a href="https://www.dgl.ai/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>dgllife</td>
            <td>0.3.2</td>
            <td class="text-center"><a href="https://lifesci.dgl.ai/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>pytorch</td>
            <td>2.1.1</td>
            <td class="text-center"><a href="https://pytorch.org/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>torchvision</td>
            <td>0.16.1</td>
            <td class="text-center"><a href="https://pytorch.org/vision/stable/index.html" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a>
            </td>
        </tr>
        <tr>
            <td>reportlab</td>
            <td>4.0.8</td>
            <td class="text-center"><a href="https://www.reportlab.com/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>scikit-learn</td>
            <td>1.3.2</td>
            <td class="text-center"><a href="https://scikit-learn.org/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>pandas</td>
            <td>2.1.4</td>
            <td class="text-center"><a href="https://pandas.pydata.org/" class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        <tr>
            <td>chemprop</td>
            <td>1.6.1</td>
            <td class="text-center"><a href="https://github.com/chemprop/chemprop"
                   class="source-link">
                <button type="button" class="btn btn-sm btn-outline-primary">Explore</button>
            </a></td>
        </tr>
        </tbody>
    </table>
</div>

## References

1. Yang,K., Swanson,K., Jin,W., Coley,C., Eiden,P., Gao,H., Guzman-Perez,A., Hopper,T., Kelley,B., Mathea,M., et al.
   (2019) Analyzing Learned Molecular Representations for Property Prediction. J. Chem. Inf. Model., 59,
   3370–3388. https://doi.org/10.1021/acs.jcim.9b00237

2. Heid,E., Greenman,K.P., Chung,Y., Li,S.-C., Graff,D.E., Vermeire,F.H., Wu,H., Green,W.H. and McGill,C.J. (2023)
   Chemprop: Machine Learning Package for Chemical Property Prediction
   Chemistry. https://doi.org/10.26434/chemrxiv-2023-3zcfl

3. Begoli,E., Bhattacharya,T. and Kusnezov,D. (2019) The need for uncertainty quantification in machine-assisted
   medical decision making. Nat. Mach. Intell., 1, 20–23. https://doi.org/10.1038/s42256-018-0004-1

4. Amini,A., Schwarting,W., Soleimany,A. and Rus,D. (2020) Deep evidential regression. Adv. Neural Inf. Process. Syst.,
   33, 14927–14937.

5. Gal,Y. and Ghahramani,Z. (2016) Dropout as a Bayesian Approximation: Representing Model Uncertainty in Deep Learning.
   In Proceedings of The 33rd International Conference on Machine Learning. PMLR, pp. 1050–1059.

6. Wen,M. and Tadmor,E.B. (2020) Uncertainty quantification in molecular simulations with dropout neural network
   potentials. Npj Comput. Mater., 6, 1–10. https://doi.org/10.1038/s41524-020-00390-8