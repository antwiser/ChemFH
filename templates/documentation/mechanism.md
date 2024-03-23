# Frequent Hitter Mechanisms

---

## <i class="bi bi-caret-right-fill"></i> Colloidal aggregators

Colloidal aggregation is a major cause of false positives observed in High Throughput Screening (HTS) across various
assay formats(1). Unlike micelle formation, the critical aggregation concentration (CAC) governs the assembly of these
aggregators. Compound aggregation typically initiates above the CAC concentration, forming aggregators with a radius of
approximately 30–600 nm(2). These colloidal aggregators can nonspecifically bind to protein surfaces, inducing local
protein unfolding and often destabilizing or denaturing enzymes (Figure 1)(3). Proteins with lower stability are more
susceptible to inhibition by colloidal aggregates(4). The binding affinity between aggregators and proteins is strong,
with corresponding KD values in the picomolar range. Larger aggregators with radii >250 nm often exhibit poorer
absorption. In the diffusion of colloid-forming dye Evans Blue, molecules can easily pass through membranes in their
free state but face difficulty in the aggregated state(5). In addition to large aggregators, certain compounds can form
smaller 'nonclassical' aggregators capable of modulating protein functions(6). Some drugs, like felodipine and benzyl
benzoate, can form aggregators in enzymatic assays at micromolar concentrations. However, these drugs are typically
optimized for potency and selectivity and are used at nanomolar concentrations, well below their CACs. Therefore, the
common strategy in early drug discovery is to filter out colloidal aggregators as early as possible.

**Results interpretation:** Category 0: Non-aggregators; Category 1: aggregators. The output value is the probability of
being aggregators, within the range of 0 to 1.

![](../static/explanation/img/aggregator.png ':class=content-img')

<center><strong>Figure 1.</strong> Colloidal aggregators: when their concentration is higher than their critical aggregation concentration (CAC), compounds start forming aggregators. Colloidal aggregators nonspecifically bind proteins, thus leading to protein destabilization or denaturation.</center>

## <i class="bi bi-caret-right-fill"></i> FLuc inhibitors

Luciferase-based bioluminescence assays are widely applied in HTS because of their high sensitivity(7). Firefly
luciferase (FLuc), obtained from the firefly Photinus pyralis, is one of the most frequently used luciferase types in
HTS. A complicated enzymatic mechanism is responsible for the light reaction of FLuc, where the substrates ATP and
D-luciferin are used by FLuc to form a luciferyl-adenylate intermediate. This intermediate then receives nucleophilic
attack by molecular oxygen, and the subsequent displacement of AMP results in the formation of an unstable dioxetanone
compound. Finally, dioxetanone breaks down spontaneously to form CO2 and oxylucifein, which results in the emission of
photons, generating both green and red light(8). Given its unique catalysis mechanism, FLuc is used in a variety of HTS
bioluminescence assays, especially in assays used to study gene expression at the transcriptional level(9). However, the
inhibitions of FLuc by unexpected FLuc inhibitors produce interference to HTS assays. There are two major types of FLuc
interferences in fluorescent detection: specific inhibition (direct inhibition of the enzyme, Figure. 2) and nonspecific
interference (enzyme denaturation or attenuation of the luminescent signal through photonic processes, such as
absorbance or the inner-filter effect).

**Results interpretation:** Category 0: Non-FLuc inhibitors; Category 1: FLuc inhibitors. The output value is the
probability of being FLuc inhibitors, within the range of 0 to 1.

![](../static/explanation/img/luciferae%20inhibitor.png ':class=content-img')

<center><strong>Figure 2.</strong> Fluorescent compounds: the active compound binds with the targets that absorb and emit light, 
whereas the fluorescent compound can exhibit similar absorption and emission spectra.</center>

## <i class="bi bi-caret-right-fill"></i> Blue/Green fluorescence

Fluorescence is the process by which a molecule, called a fluorophore or fluorescent dye, absorbs a photon of light,
exciting an electron to a higher energy state(10). Fluorophores have many applications, such as enzyme substrates,
labels for biomolecules, cellular stains, and environmental indicators(11). However, the appearance of fluorescent
compound would produce interference to related HTS assays. There are two main mechanisms responsible for fluorescent
assay interferences: quenching and autofluorescence. Quenching refers to any process that attenuates the intensity of
emitted fluorescent light. Molecular quenching is a common type of interference observed in fluorescence and luciferase
reporter enzyme assays. Autofluorescence is the natural emission of light by a given substance. Compounds whose
fluorescence region overlaps with the region of detection fluorophore are more likely to cause interferences (Figure 3)(
10).

**Results interpretation:** Category 0: Non-Blue/Green fluorescence; Category 1: Blue/Green fluorescence. The output
value is the probability of being fluorescence, within the range of 0 to 1.

![](../static/explanation/img/Fluorescent.png ':class=content-img')

<center><strong>Figure 3.</strong> Fluorescent compounds: the active compound binds with the targets that absorb and emit light, whereas the fluorescent compound can exhibit similar absorption and emission spectra.</center>

## <i class="bi bi-caret-right-fill"></i> Reactive compounds

Chemical reactive compounds (CRCs) typically result in the chemical modification of reactive protein residues or, less
frequently, the modification of nucleophilic assay reagents. The mechanisms related to chemical reactivity interference
are complicated. Moreover, some compounds are reactively inactive against targets but could be converted into CRCs in
cells through cellular metabolic processes (12). Therefore, the exploration of CRCs is challenging. Two common typical
interference mechanisms of CRCs include the redox cycling generation of hydrogen peroxide (H2O2) and electrophilic
functionality. Compounds that generate H2O2 by redox cycling in reducing agents can modulate the activity of their
targets by the oxidization of cysteine residues and, therefore, would be categorized as false-positives(13–15) . For
example, caspase-3 is an attractive therapeutic target for the treatment of diseases involving dysregulated apoptosis.
The mechanism of caspase-3 inactivation by isoquinoline-1,3,4-trione derivatives is the result of reactive oxygen
species, where the binding site of the catalytic cysteine is oxidized to sulfonic acid (–SO3H) (Figure 4)(16).

**Results interpretation:** Category 0: Non-reactive; Category 1: reactive. The output value is the probability of being
reactive compounds, within the range of 0 to 1.

![](../static/explanation/img/reactive.png ':class=content-img')

<center><strong>Figure 4.</strong> Chemical reactivity: an inhibitor binds with the target to achieve target inhibition, whereas the chemical reactive compound generates H2O2, which can oxidize cysteine residues, thus nonspecifically and promiscuously inhibiting target activity.</center>

## <i class="bi bi-caret-right-fill"></i> Promiscuous compounds

The promiscuity of compounds is a broad concept. Here, promiscuous compounds refer to those binding specifically to
different macromolecular targets. Unlikely for a molecule to occupy multiple, diverse binding sites through simple,
noncovalent interactions typical in specific screening hits, a promiscuous activity profile usually implies a
nonspecific mode of action(17). These interactions may involve unintended targets, triggering adverse reactions and
safety issues(18).

**Results interpretation:** Category 0: Non-promiscuous; Category 1: promiscuous. The output value is the probability of
being promiscuous, within the range of 0 to 1.

![](../static/explanation/img/Promiscuous.png ':class=content-img')

<center><strong>Figure 5.</strong> Promiscuous Compounds: Compounds that Specifically Bind to Different Macromolecular Targets.</center>

## <i class="bi bi-caret-right-fill"></i> Other assay interferences

Homogeneous proximity assays are widely employed in high-throughput screening of small-molecule libraries for drug and
probe discovery(19). Representative technologies include the amplified luminescent proximity homogeneous assay (ALPHA,
trademarked by PerkinElmer and informally referred to as Alpha), Förster/fluorescence resonance energy transfer (FRET),
time-resolved FRET (TR-FRET) (Figure 6). While these detection technologies are highly valuable, they are susceptible to
various technology-related compound-mediated interferences, particularly signal attenuation (e.g., through quenching,
inner-filter effects, light scattering), signal emission (e.g., autofluorescence), and disruption of affinity capture
components (such as affinity tags and antibodies), leading to false-positive results in the assays(20, 21).

**Results interpretation:** Category 0: Non-assay interferences; Category 1: assay interferences. The output value is
the probability of being assay interferences, within the range of 0 to 1.

![](../static/explanation/img/Other.png ':class=content-img')

<center><strong>Figure 6.</strong> Illustration of Interference in Förster/Fluorescence Resonance Energy Transfer 
(FRET) Method.</center>

## References

1. Auld,D.S., Inglese,J. and Dahlin,J.L. (2004) Assay Interference by Aggregation. In Markossian,S., Grossman,A.,
   Brimacombe,K., Arkin,M., Auld,D., Austin,C., Baell,J., Chung,T.D.Y., Coussens,N.P., Dahlin,J.L., et al. (eds), Assay
   Guidance Manual. Eli Lilly & Company and the National Center for Advancing Translational Sciences, Bethesda (
   MD). http://www.ncbi.nlm.nih.gov/pubmed/28749639

2. Duan,D., Torosyan,H., Elnatan,D., McLaughlin,C.K., Logie,J., Shoichet,M.S., Agard,D.A. and Shoichet,B.K. (2017)
   Internal Structure and Preferential Protein Binding of Colloidal Aggregates. ACS Chem. Biol., 12,
   282–290. https://doi.org/10.1021/acschembio.6b00791

3. Coan,K.E.D. and Shoichet,B.K. (2008) Stoichiometry and Physical Chemistry of Promiscuous Aggregate-Based
   Inhibitors. J. Am. Chem. Soc., 130, 9606–9612. https://doi.org/10.1021/ja802977h

4. Torosyan,H. and Shoichet,B.K. (2019) Protein Stability Effects in Aggregate-Based Enzyme Inhibition. J. Med. Chem.
   , 62, 9593–9599. https://doi.org/10.1021/acs.jmedchem.9b01019

5. Owen,S.C., Doak,A.K., Ganesh,A.N., Nedyalkova,L., McLaughlin,C.K., Shoichet,B.K. and Shoichet,M.S. (2014)
   Colloidal Drug Formulations Can Explain “Bell-Shaped” Concentration–Response Curves. ACS Chem. Biol., 9,
   777–784. https://doi.org/10.1021/cb4007584

6. Structural Basis of Small-Molecule Aggregate Induced Inhibition of a Protein–Protein Interaction | Journal of
   Medicinal Chemistry.

7. Auld,D.S. and Inglese,J. (2004) Interferences with Luciferase Reporter Enzymes. In Markossian,S., Grossman,A.,
   Brimacombe,K., Arkin,M., Auld,D., Austin,C., Baell,J., Chung,T.D.Y., Coussens,N.P., Dahlin,J.L., et al. (eds), Assay
   Guidance Manual. Eli Lilly & Company and the National Center for Advancing Translational Sciences, Bethesda (
   MD). http://www.ncbi.nlm.nih.gov/pubmed/27478246

8. Thorne,N., Inglese,J. and Auld,D.S. (2010) Illuminating Insights into Firefly Luciferase and Other Bioluminescent
   Reporters Used in Chemical Biology. Chem. Biol., 17, 646–657. https://doi.org/10.1016/j.chembiol.2010.05.012

9. Auld,D.S., Southall,N.T., Jadhav,A., Johnson,R.L., Diller,D.J., Simeonov,A., Austin,C.P. and Inglese,J. (2008)
   Characterization of Chemical Libraries for Luciferase Inhibitory Activity. J. Med. Chem., 51,
   2372–2386. https://doi.org/10.1021/jm701302v

10. Simeonov,A. and Davis,M.I. (2004) Interference with Fluorescence and Absorbance. In Markossian,S., Grossman,A.,
    Brimacombe,K., Arkin,M., Auld,D., Austin,C., Baell,J., Chung,T.D.Y., Coussens,N.P., Dahlin,J.L., et al. (eds), Assay
    Guidance Manual. Eli Lilly & Company and the National Center for Advancing Translational Sciences, Bethesda (
    MD). http://www.ncbi.nlm.nih.gov/pubmed/26844333

11. Lavis,L.D. and Raines,R.T. (2008) Bright Ideas for Chemical Biology. ACS Chem. Biol., 3,
    142–155. https://doi.org/10.1021/cb700248m

12. Yang,Z.-Y., He,J.-H., Lu,A.-P., Hou,T.-J. and Cao,D.-S. (2020) Frequent hitters: nuisance artifacts in
    high-throughput screening. Drug Discov. Today, 25, 657–667. https://doi.org/10.1016/j.drudis.2020.01.014

13. Smith,G.K., Barrett,D.G., Blackburn,K., Cory,M., Dallas,W.S., Davis,R., Hassler,D., McConnell,R., Moyer,M. and
    Weaver,K. (2002) Expression, Preparation, and High-Throughput Screening of Caspase-8: Discovery of Redox-Based and
    Steroid Diacid Inhibition. Arch. Biochem. Biophys., 399, 195–205. https://doi.org/10.1006/abbi.2002.2757

14. Bova,M.P., Mattson,M.N., Vasile,S., Tam,D., Holsinger,L., Bremer,M., Hui,T., McMahon,G., Rice,A. and Fukuto,J.M.
    (2004) The oxidative mechanism of action of ortho-quinone inhibitors of protein-tyrosine phosphatase α is mediated
    by hydrogen peroxide. Arch. Biochem. Biophys., 429, 30–41. https://doi.org/10.1016/j.abb.2004.05.010

15. Redox Regulation of Cdc25B by Cell-Active Quinolinediones | Molecular Pharmacology.

16. Du,J.-Q., Wu,J., Zhang,H.-J., Zhang,Y.-H., Qiu,B.-Y., Wu,F., Chen,Y.-H., Li,J.-Y., Nan,F.-J., Ding,J.-P., et al.
    (2008) Isoquinoline-1,3,4-trione Derivatives Inactivate Caspase-3 by Generation of Reactive Oxygen Species*. J.
    Biol. Chem., 283, 30205–30215. https://doi.org/10.1074/jbc.M803347200

17. Bruns,R.F. and Watson,I.A. (2012) Rules for Identifying Potentially Reactive or Promiscuous Compounds. J. Med.
    Chem., 55, 9763–9772. https://doi.org/10.1021/jm301008n

18. Senger,M.R., Fraga,C.A.M., Dantas,R.F. and Silva,F.P. (2016) Filtering promiscuous compounds in early drug
    discovery: is it a good idea? Drug Discov. Today, 21, 868–872. https://doi.org/10.1016/j.drudis.2016.02.004

19. Coussens,N.P., Auld,D., Roby,P., Walsh,J., Baell,J.B., Kales,S., Hadian,K. and Dahlin,J.L. (2004)
    Compound-Mediated Assay Interferences in Homogeneous Proximity Assays. In Markossian,S., Grossman,A., Brimacombe,K.,
    Arkin,M., Auld,D., Austin,C., Baell,J., Chung,T.D.Y., Coussens,N.P., Dahlin,J.L., et al. (eds), Assay Guidance
    Manual. Eli Lilly & Company and the National Center for Advancing Translational Sciences, Bethesda (
    MD). http://www.ncbi.nlm.nih.gov/pubmed/32045180

20. Aldrich,C., Bertozzi,C., Georg,G.I., Kiessling,L., Lindsley,C., Liotta,D., Merz,K.M., Schepartz,A. and Wang,S.
    (2017) The Ecstasy and Agony of Assay Interference Compounds. J. Med. Chem., 60,
    2165–2168. https://doi.org/10.1021/acs.jmedchem.7b00229

21. Dahlin,J.L. and Walters,M.A. (2016) How to Triage PAINS-Full Research. Assay Drug Dev. Technol., 14,
    168–174. https://doi.org/10.1089/adt.2015.674