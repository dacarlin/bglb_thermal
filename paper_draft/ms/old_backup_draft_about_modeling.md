# Experimental measurements of protein expression, thermal stability, and enzyme kinetic constants for 128 mutants of a glycoside hydrolase enable discovery of structural features correlated to protein stability

## Author contributions (alphabetical by last name)

+ Bowen Bethards [1]: designed mutants, characterized mutants
+ Dylan Alexander Carlin [2]: designed mutants, cloning, designed experiments, wrote software used in analysis, analyzed data, modeling, machine learning, wrote paper
+ Ryan Caster [1]: characterized mutants
+ Bill Chan [1]: designed mutants, characterized mutants, analyzed data, contributed to paper
+ Natalie Damrau [1]: designed mutants, characterized mutants
+ Siena Hapig-Ward [1]: designed mutants, characterized mutants, analyzed data, contributed to paper
+ Mary Riley [1]: designed mutants, characterized mutants
+ Justin B. Siegel [1,3,4]: PI, provided lab, materials, designed experiments, paper

## Author affiliations:

1. Genome Center, University of California, Davis CA, USA
2. Biophysics Graduate Group, University of California, Davis CA, USA
3. Department of Chemistry, University of California, Davis CA, USA
4. Department of Biochemistry & Molecular Medicine, University of California, Davis CA, USA

**Subject areas**: biophysics, computational biology, protein engineering

**Keywords**: enzyme, Rosetta, thermal stability, machine learning

---

## Abstract

The importance of enzymes in biotechnology makes the accurate modeling of enzyme thermal stability and soluble expression in heterologous hosts such as *Escherichia coli* an important goal for the protein modeling community. Previous efforts at accurate computational prediction, often using have been partially successful, but are hampered by the lack of large data sets on which algorithms can be trained. Here, we report soluble protein expression, thermal stability (Tm), and Michaelis-Menten constants kcat and KM for 128 individually-purified variants of a family 1 glycoside hydrolase enzyme. We demonstrate the utility of this data set by training a machine learning algorithm to predict soluble protein expression and Tm from readily-calculable features. With the experimental data set and analyses carried out in this study, we discover unexpected functional effects of mutations that shed light on enzyme function, as well as provide a critical assessment of current predictive algorithms and goals for improvement.

---

## Introduction

[Option 1: Human health pitch]

The ability to predict the stability effects of mutations to enzymes is a major goal of the protein modeling community, with dozens of software packages available [molec_mech]. The major goal of these efforts has been to determined the stability changes induced by single point mutations, as mutations that lead to protein misfolding and misfunction cause disease in humans [molec_mech].

However, accurate prediction of enzyme functional envelope from readily-calculable structural features remains an unmet challenge.

[Human health pitch]

[Option 2: Biotech pitch]

The accurate modeling of enzyme thermal stability and soluble expression in heterologous hosts such as *Escherichia coli* an important goal for the protein modeling community, with many tools available for making predictions from sequence, structure, or a combination of sequence, structure, and experimental data [molec_mech]. Yet, despite the importance of these goals, the prediction of the stability effect of point mutations in proteins has not been fully achieved, and neither has the prediction of the expression of sequences in heterologous hosts.  

[End biotech pitch]

[Option 3: function-stability tradeoff pitch]

The idea that functional proteins, such as enzymes, must trade stability in order to carry out their functions is widely believed to be a general principle of enzymatic catalysis [Beadle, Tokuriki, RNase paper]. In particular, enzyme active sites, such as the one in BglB, often position charged residues in close proximity in order to present optimal geometry for catalysis. [Fersht] In an enzyme active site such as those found in family 1 glycoside hydrolases, two glutamate residues, one functioning as a nucleophile and the other as an acid/base in a classical Kosland double-displacement mechanism, are positioned at 5.5 A.

In order to test this hypothesis, we designed 128 point mutants of BglB and individually characterized the soluble protein expression, thermal stability, Kcat, KM, and kcat/KM of each mutant. We present this data set, along with our analysis, here.

[Beadle] Structural Bases of Stability–function Tradeoffs in Enzymes
[Tokuriki] How Protein Stability and New Functions Trade Off
[RNase paper]

Dozens of software packages are available for prediting the stability effects of mutations [molec_mech]. Some use features derived from evolutionary information and sequence alone [SIFT], while some incorporate structural information [PolyPhen, AUTO-MUTE, PopMuSiC]. Some, like FoldX, use force-field–based molecular modeling combined with a statistical potential that has been trained on experimental data for point mutations [FoldX]. Other previous work has incorporated machine learning algorithms [Masso, Berliner], using the metrics generated by the various software packages as features. None, however, has been independently evaluated using an independent, curated training set.  

Data sets that have been used to train the software packages above have three primary weaknesses. Some are large are constructed without regard for standardization [ProTherm], reducing the quantitative accuracy and usefulness of predictions made by these methods. The largest data sets for which melting temperature $Tm$ is explicitly measured have been around 30 mutants [cite]. Other studies have collected large amounts of data, but it suffers from convolution of separate/othoganal [we will show that kcat and tm are not correlated] parameters/fuzzy measuremetns [cite: Romero and others]. The result is that existing data sets either contain specific measurements for a small number of enzymes, or fuzzy measurements of large numbers of variants. This is likely due to the immense cost of generating and characterizing large numbers of variants under standardized conditions.

A further problem is that high-quality data sets that are quantitative and large have already been used to refine protein modeling algorithms and thus in some sense have been seen by the algorithms before. Thus, the accurate modeling of stability effects of point mutations, while an important goal, remains a substantial challenge for protein stability prediction algorithms.

No study has yet combined both a standardized approach to experimental characterization, and produced enough mutants to allow a sufficiently large training set to evaluate predictive ability. We have previously reported on a data set containing of Michaelis-Menten constants kcat, KM, and kcat/KM for 100 mutants of a glycoside hydrolase. Using a constrained statistical learning technique (elastic net), we evaluated 50 structural features derived from modeling using RosettaDesign [Carlin]. Building on our previous work, we have experimentally determined kinetic constants for a further 28 mutants, as well as soluble protein expression and thermal stability for the entire data set of 128 mutants.

Here, we show that this data set can be used as training data to evaluate the utility of computational features generated by protein modeling software Rosetta and FoldX in predicting the soluble expression of the 128 mutants in our data set. We then show the use of the experimental data to build a predictive model of soluble enzyme expression and melting temperature using machine learning algorithms. We show that, while models for predicting melting temperature require larger data sets to become feasible, we are able to apply our model to a blind test set of 10 mutants, on which we achieved a predictive accuracy of X percent, representing a significant improvement of force-field based approaches alone and a path forward for the field.

![](figure_1.png)

**Figure 1: Overview of BglB structure showing positions mutated in this study**: PyMOL rendering of modeled BglB in complex with pNPG (top panel) showing the Cɑ of the 60 sequence positions selected for this study (purple) and the modeled transition-state structure (green). Reaction scheme of the hydrolysis of pNPG by BglB (bottom panel).

## Methods

### Mutant selection by molecular modeling in Foldit

A crystal structure of recombinant BglB in complex with the substrate analog 2-deoxy-2-fluoro-α-D-glucopyranose (PDB ID: 2JIE) was used to identify the substrate binding pocket and the catalytic residues. A macromolecular model of BglB was created using Rosetta. Functional constraints were used to define catalytic distances, angles, and dihedrals among 4-nitrophenyl-β-D-glucoside, E164, E353, and Y295. The structure was then loaded into Foldit, a graphical user interface to Rosetta. Point mutations to the protein were modeled, and a subset were chosen by students learning about molecular modeling. Generally, the designs had energies no more than 5 Rosetta energy units higher than the native structure, but no other design constraints were used to select mutations.

### Molecular cloning and mutagenesis

A gene coding for BglB was codon-optimized for *Escherichia coli* and synthesized as linear dsDNA (Life Technologies). The construct was cloned into a pET29b+ vector linearized with XhoI and NdeI (NEB) using Gibson assembly master mix (NEB) and verified by Sanger sequencing (Operon). Kunkel mutagenesis was used to generate designed variants of BglB via an automated cloud laboratory (Transcriptic), which were then verified by Sanger sequencing (Operon, GenScript).

### Protein production and purification

Sequence-verified constructs were transformed into chemically competent *Escherichia coli* BLR. Single colonies were used to inoculate 5 mL Terrific Broth (TB) (Fisher) in 50 mL Falcon tubes (Fisher) with breathable seals (Fisher). After incubation at 37 C with shaking at 300 RPM for 24 hours, cells were pelleted and media replaced with 5 mL TB containing 1 mM isopropyl β-D-1-thiogalactopyranoside (IPTG). After incubation at 18 C with shaking at 300 RPM for 24 hours, cells were pelleted, resuspended in enzyme storage buffer (50 mM HEPES, 150 mM sodium chloride, 25 mM EDTA, pH 7.50) and lysed with BugBuster protein extraction reagent (Novogen).

His-tagged BglB proteins were purified via immobilized metal ion affinity chromatography using 50 µL bed volume of Ni-NTA resin (Thermo) and eluted in 300 µL enzyme storage buffer (wash buffer was the same as enzyme storage buffer except omitting EDTA). Protein expression was assessed using 4-12% gradient SDS-PAGE (Life Technologies).

### Enzyme assays and data analysis

For thermal stability assays, triplicate aliquots of freshly-purified proteins in enzyme storage buffer at concentrations ranging from 0.01 to 0.1 mg/mL were incubated at 8 temperatures in 2.5 C increments between 30 C and 50 C in a thermal cycler (BioRad). After 30 minutes, proteins were transferred to assay plates containing 100 mM pNPG in enzyme storage buffer. The rate of pNP production was monitored over 60 minutes and a linear fit to the data was determined by Gen5 software.

For each BglB variant, triplicate rates at 8 temperatures were normalized to the [0,1] interval and fit to the logistic equation `1/(1+exp(-k*(x-x0)))` using least squares (accessed via SciPy's `curve_fit`) to determine the parameters Tm (midpoint) and k (kurtosis) and one-standard-deviation errors from fitted values.

Michaelis-Menten parameters for each mutant were determined as described previously [Carlin]. Mutants for whom Michaelis-Menten parameters were previously reported use values from the publication.

## Computational modeling and feature generation using Rosetta and FoldX

Three separate approaches to modeling the 128 mutants were taken in this study. In the first, 100 molecular models of each mutant enzyme were generated using the Rosetta Molecular Modeling Suite by Monte Carlo optimization of total system energy using enzyme design protocols as reported in our previous work [Carlin]. The lowest 10 of total system energy were selected for feature generation and a total of 50 features were calculated for each protein model. In a second approach, mutants were generated and scored using benchmarked [Kellogg] Rosetta ddG application, which generates 50 models and averages the best 3 models to provide a feature set of 13 features. In a third approach, models were generated and scored using FoldX PSSM application, which generates a feature set of 13 features per model. The lowest energy model was selected for feature generation.

### Machine learning methods

Elastic net for Tm

Random forests for expression

# Results

A total of 128 proteins, not including batched replicates of wild type BglB produced with each batch, were produced, purified, and characterized in this study. Of 128 proteins synthesized and assayed, 86 (67%) express and purify as soluble protein. The remaining 42 mutants did not appear on a gel after at least 2 independent production attempts.

Of the 86 solubly-expressed mutants, melting temperature for 67 was determined by fitting rate data collected at 8 temperatures from 30 to 50 C to the logistic equation. For the the remaining 21 mutants that expressed as soluble protein,  could not be used to determine Tm because the kinetic constants for these mutants are below our limit of detection.

The wild type BglB replicates had an average melting temperature of 39.6 C across assay data from 8 biological replicates. The average melting temperature was 39.0±1.7 C, with a range of 32.9 C to 45.3 C. Of 62 mutants for which Tm was determined, Tm for 37 of them falls within 1 degree C of the wild type Tm (58%). Of the remaining 49 Tm values, 18 exhibited a lower melting temperature than that of native BglB and 7 displayed a higher melting temperature. The hottest Tm observed in this data set is N404C, which increased the Tm of BglB to 45.3 C (dTm=+5.4 ˚C). The mutation that decreased the thermal stability of BglB the most was H178A (dTm=-7.0 ˚C).

In addition to expression and thermal stability data, we also report kinetic constants for 12 mutants whose kinetic constants had not previously been determined. Together, data for soluble expression, kinetic constants (kcat, km, kcat/km), and melting temperature are reported for 128 mutants of BglB (Supporting Information). **Figure 2** depicts the data set as a heat map, with values relative to native BglB. Based on the maximum concentration of enzyme used in our assays and colorimetric absorbance changes at the highest substrate concentration used, we estimate our limit of detection for kcat/KM to be 10 M<sup>-1</sup>min<sup>-1</sup>.

![](figure_2.png)

**Figure 2: heatmap showing the expression, Tm, kcat, KM, kcat/KM, and percent conservation for 128 mutants**. For expression (first column), a black box indicates soluble expression, and a white box indicates no expression for this mutant. The last column indicates the degree of conservation for the BglB native residue based on an alignemnt of 1,554 sequences from Pfam. For all other values, the heatmap depicts each value relative to the wild type value. For Tm, a linear scale from -6 (purple) to 6 (yellow) is used to depict deltaTm compared to wild type. For kcat, KM, and kcat/km, a diverging colormap is used, with blue values indicating lower values, and yellow indicating higher values, as indicated by the color legend. For KM, the value 1/KM is plotted, so that a higher plotted value indicates greater affinity between BglB and pNPG for the mutation.

## Sequence-stability relationships in BglB

In agreement with previous studies, mutants that did not express followed rules broadly consistent with previous work and sequence conservation, such as the large destabilizing effect of any substitution at W407. Other mutants that did not express mostly followed well-established observations of destabilizing effects, such as the introduction of proline into an alpha helix (Y166P, Q19P), the mutation of topology-defining/helix-ending proline residues (P329N), mutations to and from glycine (), the introduction of charged residues into the hydrophobic core (A236E, F72H), and extreme small-to-large mutations (A227W).

It is interesting to note that all the mutations we made in our study to the catalyic nucelophile (E353) and acid-base (E164) polar residues resulted in soluble protein, consistent with the idea that enzymes must trade structural stability for function in placing these two negatively-charged groups less than 5 Å apart (the sole exception, E164G, could be the result of an altered folding trajectory due to the conformational freedom of glycine). Notable is that all the catalytic knockouts, E164A, E164R, E353A, Y295A, Y295G, except for E164G, where a glycine is inserted in an alpha helix, result in protein that is solubly expressed in our system, supporting the belief that BglB must compromise structural stability to perform its function.

A novel finding was a nearly six degree increase in melting temperature by single point mutant N404C. The BglB crystal structure reveals a weak hydrogen bond between N404 and the backbone at L402. Molecular modeling of N404C predicts the loss of this hydrogen bond to the protein's alpha helix, in which the protein is allowed to repack into more energetically favorable states.

Similarly, the point mutation W120F resulted in a delta Tm of +1.6 C. The BglB crystal structure indicates a weak hydrogen bond between W120 and the backbone of N163, which could be construed as an unsatisfied polar interaction. The mutation to Phe maintains the structural integrity at the mutation site as well as removes the unsatisfied hydrogen bond donor to the neighboring alpha helix. The increased stability is then due to destabilization of the unfolded state, which exposes a hydrophobic Phe to bulk solvent. It is also worth noting that the mutant W120A results in no soluble protein production after 3 independent attempts, indicating that W120 plays a key role in stabilizing the protein. Previous studies have shown a similar increase in stability upon W to Y point mutations [Fulton]. Analysis of a multiple sequence alignment of 1,554 proteins from the Pfam database reveals approximately equal probability of finding Trp, Phe, or Tyr here, and less than 1% representation of all of the other amino acids combined. Thus, the experimental measurements and the evolutionary record agree that W120 plays a key role in stabilizing BglB. No major structural changes are predicted via Rosetta modeling.  

Additionally, we have additional information about R240. The mutation R240E, which plates a glutamate near the existing glumatate at E222 that forms two hydrogen bonds to R240, placing two negativly-charged residues in close proximity, results in no soluble protein under our experimental conditions after at least 2 repeated attempts to produce and purify the protein, providing evidence for the belief that the mutation R240E is destabilizing. Whereas Rosetta modeling predicts no change in deltaG (within 1 REU) for this mutation, and example of a mutation that would have been missed using Rosetta. However, FoldX assigns a ddG of 3.78, highlighting the importance of a diverse set of estimates chosen by machine learning on experimental data to train algorithms. (Mutation R240D, which places a carboxylate only 1.54 A further away, results in a soluble protein.)

Point mutation E222H had a melting temperature of 34.7 C, a nearly five degree decrease than that of native BglB. Previous studies show strong hydrogen bond interaction, 2.6 and 3.1 Å, between E222 and its neighboring R240 residue [cite]. The introduction of histidine at the mutation site causes the loss of these strong hydrogen bonds as well as the creation of electrostatic repulsion between the partially positively charged and positively charged amino acids. The cumulative effect of this mutation results in the protein's decreased stability at lower temperatures.

![](figure_3.png)

**Figure 3: Single point mutants with extreme effects on thermal stability**:  Four mutant panels are shown, sorted from left to right by delta-Tm. For each mutant, PyMOL images depict the local area of the mutation in the BglB WT protein (top) and Rosetta model of mutatation (bottom). In the bottom most panel, the assay for the mutant (black dots) is shown along with the fitted Tm (black dashed line) and the melting curve for WT BglB (green) as measured over 8 replicates in this study (the Y axis is in arbitrary units and the X axis covers the range 30-50 C).

## Statistical analysis of single feature correlation for Rosetta and FoldX

So, given our large data set, we have a set of questions.

The first question is, how good are Rosetta and FoldX are at predicting Tm? To answer this question, we modeled each of 128 mutations using accepted, validated computational protocols. We then asked, for each feature and for each data set, what is the correlation between the experimental Tm and the computational feature (as assessed by Pearson correlation coefficient).

[See Jupter notebook]

**Table 1. Correlation between Tm and structural features generated by RosettaDesign, Rosetta ddG, Rosetta FastRelax, and FoldX PSSM.** Each correlation is given as the Pearson Correlation coefficient.  

[See Jupter notebook]

**Figure X. Scatter plots depicting the relationship between Tm and each of structural features generated by RosettaDesign, Rosetta ddG, Rosetta FastRelax, and FoldX PSSM.** The experimental Tm is on the Y axis and the calculated structural feature is on the X axis.   

After performing this statistical analysis in the provided Jupyter notebook, I found that the best correlation between Tm and the features in our data set had a PCC of 0.3, a value that falls far short of our goal of accurate prediction of Tm for this data set. The results are summarized in a table below.

| Feature set | Number of features | Most-correlated feature to T<sub>m</sub> | PCC |
| ----------- | ------------------ | ---------------------------- |
| RosettaDesign | 60 | SR_2_pstat_pm | 0.31 |
| Rosetta ddG | 15 | hbond_bb_sc | 0.36 |
| Rosetta FastRelax | 20 | omega | 0.20 |
| FoldX | 23 | Backbone Hbond | 0.35 |

**Table 2: Summary of single feature correlations with protein melting temperature.** The best-correlated feature (as assessed by PCC) from each of 4 feature sets used to predict protein melting temperature.

The next question is, given that PCC values of 0.36, such as we are discovering to be the case in our data set, are not considered predictive, we think that neither Rosetta nor FoldX alone is sufficient.

## Machine learning discovers informative features correlated with Tm

What would happen if we used machine learning to build predictive models of Tm using the data set as training data? Using a cross-validation procedure, we can estimate the generalization error of our predictions. In our previous work, we used an elastic net algorithm to simultaneously select features and build a predictive model for enzyme kinetic parameters. Can the same be applied here?

# Discussion

The Rosetta Molecular Modeling Suite has been successfully used to guide the engineering of a wide range of enzyme functions. However, there has been a limited ability to benchmark its predictive power for enzyme reengineering due to the lack of a large, quantitative dataset correlating the effects of mutations to changes in kinetic parameters and protein stability. Here, we expand on a previous data set of 100 mutants and report statistical evaluation of our ability to predict the functional effects of enzyme mutations on enzyme kinetic parameters, protein expression, and thermal stability.

The data set that forms the basis of this study uncovered unexpected structure-function relationships in BglB, including mutants with significant effects on thermal stability not predicted by computational modeling and having no evolutionary precedent in available databases.

Traditionally, studies of this type have sought to design more thermostable mutants. However, we are interested more in less thermostable mutants in the context of human disease. Many inherited disorders of metabolism are caused by random point mutations in germline protein-coding genes. Disease results when the mutant protein lacks a key catalytic residue, or when protein misfolding or instability leads to an inactive form. Since two alleles of metabolic enzymes are simultaneously expressed, a precise quantitative measurement of the kinetic activity and protein stabilty is necessary for the predicition to have any relevance.

![](figure_4.png)

**Figure 4: statistical independence of protein melting temperature (Tm) and kinetic constants.** Plots of Tm versus kinetic constant for each of kcat, KM, and kcat/km. Tm values are on a linear scale (C) and kinetic constant values are on a log scale.

The finding that protein melting temperatures and kinetic constants reported in this study were not found to be correlated. This sheds doubt on studies [Romero] that convolute these parameters into specific activity measurements, for one can never be sure the functional parameter that is tuned to result in an inactive mutant. In biotechnology, however, this finding is a boon: it means that kinetics and thermal stability, and perhaps other parameters of the functional envelope, arise from independent (if overlapping) biophysical properties, and thus can be rationally modulated independently. This neatly avoids the "two-objective optimisation problem" assumed to exist by engineers seeing to maximize two parameters (such as kcat and thermal stability) independently [http://pubs.rsc.org/en/content/articlepdf/2015/cs/c4cs00351a]. Whether this finding is relevant to BglB alone is unclear, as other studies have shown a strict linear correlation between delta-G and kcat/km for RNase [http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1219248/]. However, it will be impossible to know without large data sets such as the one preseneted here.

One major reason for the development large, self-consistent data sets such as the one presented here is that they do not exist in existing databases and thus provide unbiased results. From the point of view that they are not already in the training set, as would happen if FoldX or Rosetta scoring functions are updated without publication, which often happens. As published data sets are used to refine existing software, that software accumulates a bias toward overfitting that is recapitulated whenever a new group seeks to evaluate the predictive ability of the software on old data sets. Therefore, it is of crucial importance to generate large, self-consistent data sets such as the one presented here to accurately evaluate predictive modeling algorithms.

Since to the challenges faced in experimentally characterizing large numbers of mutants, which often carries significant financial cost, are so great, computational approaches to predict the stability of enzyme mutations have in turn relied on these far from ideal data sets as training data. Previous work to develop predictive models of enzyme function using a combination of machine-learning and features derived from force-field based computational modeling of protein structure have relied in some cases on protein-protein-interaction data for large number of individual proteins, with a small number of mutations each [PPI dataset], and in others  ProTherm, a database that seeks to collect data on enzyme mutations. ProTherm has the additional problem that no attempt at standardization has been made, and as a result the experimental values are on different scales for each mutation and certainly for each protein. Previous attempts to predict the effect of point mutations on antibodies (~100 residues) have used primary sequence and amino acid properties as features \cite{21710487}.

---

### Supporting information

  1. Table of protein expression, melting temperature, kinetic constants, and statistical analysis for each of 125 mutants (CSV)
  1. Michaelis-Menten plots for each mutant for which kinetic constants are reported (ZIP file containing images)
  1. Plot of protein melting for each mutant for which Tm is reported (ZIP file containing images)
  1. Images of gel band for each protein used in this study (ZIP file or single TIFF image)
  1. Additional information about experimental procedures (text)
  1. Additional information about machine learning (text)
  1. Rosetta and FoldX input files (ZIP containing text files)

---

## References

[antibody_corr] High-throughput measurement, correlation analysis, and machine-learning predictions for pH and thermal stabilities of Pfizer-generated antibodies

[Masso] Accurate prediction of stability changes in protein mutants by combining machine learning with structure based computational mutagenesis

[Berliner] http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4170975/

[PPI dataset] http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1003592

[molec_mech] Molecular mechanisms of disease-causing mutations JBC

## Extra stuff (crap)

(Since it is accepted that a number of human diseases are caused by defects in particular enzymes of metabolism (i.e. kcat too high, KM too low for proper function). It would be interesting to investigate whether the opposite is true: whether mutations that increase turnover number or decrease dissociation constants also result in diseases. Perhaps diseases of excess are as common as diseases where a defective enzyme is to blame.)

In any case, our previous work sampled functional space in a way that captured almost entirely constants of lower magnitude than the natural sequence chosen as our baseline (it cannot be helped to think that if we had instead begun with a mutant of approximately half the kcat of the WT, we should have been able to see many mutations that increase the kcat). However, in the set of protein melting temperatures reported here, 25 mutants have a Tm that is greater than WT (39.6 C), while 32 have a Tm that is lower (0.78:1 ratio).

The data generated here provides a standardized, quantitatively-comparable data set for the assessment of current algorithms for predicting protein stability and soluble expression in heterologous hosts, which will lead to the improvement of computational tools, with far-reaching applications in biotechnology and human health.

The key role of enzymes in performing specific chemical transformations under conditions dictated by the process environment in biotechnology makes the effective modulation of enzyme functional envelope of great importance in biotechnology. On the other side of the functional spectrum, many human diseases are caused by single point mutations to enzymes that result in nonfunctional proteins. In this case, it would be of great benefit to use cheap, fast, and plentiful genomic information to predict the functional effects of SNPs.

Together, the importance of enzymes in biotechnology and human disease make the accurate modeling of enzyme stability an important goal of the protein modeling community. Accurate prediction of a point mutation's effect on enzyme stability would unlock rational protein engineering approaches, where the information could be used immediately to rationally engineer an enzyme's functional envelope for a desired situation as has been previously explored \cite{22575958}. Furthermore, understanding the changes in enzyme stability that occur upon point mutations would provide huge insight into understanding inherited diseases of metabolism [cite], cancer [cite], as well as mechanisms by which bacteria become resistant to antibiotics, which is lately a public health menace [cite].
