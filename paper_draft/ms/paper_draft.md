# Investigation of function-stability tradeoffs in enzymes using experimental measurements of thermal stability and enzyme kinetic constants for 125 glycoside hydrolase mutants

## Author contributions (alphabetical by last name)

+ Bowen Bethards [1]: designed mutants, characterized mutants
+ Dylan Alexander Carlin [2]: designed mutants, molecular cloning, mutagenesis, designed experiments, wrote software used in analysis, analyzed data, modeling, machine learning, wrote paper
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

**Keywords**: enzyme, computational biology, biophysics

---

## Abstract

The idea that functional proteins such as enzymes trade off thermodynamic stability for functional properties is widely debated in the field of protein science. Here, we report soluble protein expression in E. coli, thermal stability, and Michaelis-Menten constants kcat, KM, and kcat/KM for 125 single point mutants of a family 1 glycoside hydrolase enzyme. Statistical analyses reveal that correlation between Tm and kcat (PCC -0.04), KM (-0.14), and kcat/KM (-0.02) in the BglB system. The data set reported here allowed us perform a general analysis of the correlation between thermal stability and kinetic constants for 125 individual single steps in sequence-functional space immediately adjacent to the native BglB sequence (Hamming distance of 1), providing a system-level look at the structure-stability-functional relationships between residues in the BglB protein, and revealing that sequence-functional tradeoffs are not a general rule, and are strongly dependent on structural features.

---

## Introduction

Whether enzymes must trade thermodynamic stability for functional properties such as pre-ordered active sites in order to carry out their functions is widely debated in the field of protein biophysics. Are stability-function tradeoffs a general principle of enzymatic catalysis [Beadle, Tokuriki, RNase paper], or are protein stability and functional parameters independent?

Enzyme active sites often position charged residues in close proximity, creating unfavorable electrostatic interactions between residues in order to present pre-ordered geometry for catalysis [Fersht]. In family 1 glycoside hydrolases, the carboxyl group of two glutamate residues, one functioning as a nucleophile and the other as an acid/base in a Koshland double-displacement mechanism, are positioned at 4.0 A in a crystal structure of a BglB-inhibitor complex. BglB relies on the proximity of the pair of glutamate residues to cyclicly preturb the pKa of E164 during catalysis, allowing it to act as a general acid in the glycosylation of the enzyme, and a general base in the product release step [Harris].

Ferst 1991 found a relationship between the thermodynamic stability of mutatns and thier kinetic activy in barnase. They found that mutation of positively-charged residues in the enzyme active site increased the stability of the protein, in one case by X kcal/mol. Similarly, Matthews' work on T4 lysozyme has shown that some evidence for a relationship between protein stability and function, providing specific activity measurements over 4 orders of magnitude for about 20 single single mutations of T4 lysozyme, as well as Tm and ddG measurements for each mutant (Shoichet_1995). However, in both these cases, only a small number of sequence positions are covered by mutagenesis experiments, in the case of Fersht 1991, only X positions, in the case of Stoichet, only 5 sequence positions covered.

In order to test the hypothesis that individual residue identities in the BglB sequence are a trade off between function and thermal stability, we designed 125 point mutants covering 68 sequence positions of BglB and individually characterized the soluble protein expression, thermal stability, kcat, KM, and kcat/KM of each mutant. Statistical analyses reveal that Tm and kinetic parameters are independent, as assessed by Pearson correlation. For kcat (PCC -0.04), KM (-0.14), and kcat/KM (-0.02), and it should be noted that all correlations have a negative algebraic sign.

![](figure_1.png)

**Figure 1: Overview of BglB structure showing positions mutated in this study**: PyMOL rendering of modeled BglB in complex with pNPG (top panel) showing the Cɑ of the 60 sequence positions selected for this study (purple) and the modeled transition-state structure (green). Reaction scheme of the hydrolysis of pNPG by BglB (bottom panel).

## Methods

### Mutant selection

A crystal structure of recombinant BglB in complex with the substrate analog 2-deoxy-2-fluoro-α-D-glucopyranose (PDB 2JIE) was used to identify the substrate binding pocket and the catalytic residues, and to build macromolecular models of BglB using Rosetta. In RosettaDesign simulations, functional constraints were used to define catalytic distances, angles, and dihedrals among a parameterized representation of 4-nitrophenyl-β-D-glucoside, and protein sidechains E164, E353, and Y295, as reported previously.

To select a subset of the 8880 possible single point mutations to the native BglB sequence that could be tested, three approaches were taken. In the first, all residues within 12 A of the modeled pNPG were individually mutated to alanine. Second, a small set of mutants were chosen at random from SNP-accessible mutations within 12 A of the active site. Third, the structure was loaded into Foldit, a graphical user interface to Rosetta, point mutations to the protein were modeled, and a subset were chosen by students learning about molecular modeling. Generally, the student-designed sequences had energies no more than 5 Rosetta energy units higher than the native structure, but no other design constraints were used to select mutations.

### Molecular cloning and mutagenesis

A gene sequence coding for BglB was codon-optimized for expression in *Escherichia coli* and synthesized as linear double-stranded DNA (Life Technologies). The construct was cloned into a pET29b(+) vector linearized with XhoI and NdeI (NEB) using Gibson assembly master mix (NEB) and verified by Sanger sequencing (Operon). Kunkel mutagenesis was used to generate designed variants of BglB via an automated mutagenesis procedure carried out by robotic automation in a remote laboratory (Transcriptic). Individual plasmid constructs were verified by Sanger sequencing (Operon, GenScript).

### Protein production and purification

Sequence-verified constructs were transformed into chemically competent *Escherichia coli* BLR. Single colonies were used to inoculate 5 mL Terrific Broth (TB) (Fisher) in 50 mL Falcon tubes (Fisher) with breathable seals (Fisher). After incubation at 37 C with shaking at 300 RPM for 24 hours, cells were pelleted and media replaced with 5 mL TB containing 1 mM isopropyl β-D-1-thiogalactopyranoside (IPTG). After incubation at 18 C with shaking at 300 RPM for 24 hours, cells were pelleted, resuspended in enzyme storage buffer (50 mM HEPES, 150 mM sodium chloride, 25 mM EDTA, pH 7.50) and lysed with BugBuster protein extraction reagent (Novogen).

After clarification of lysis mixture by centrifugation at 14,700 RPM for 30 minutes, His-tagged BglB proteins were purified via immobilized metal ion affinity chromatography using 50 µL bed volume of Ni-NTA resin (Thermo) and eluted in 300 µL enzyme storage buffer (wash buffer was the same as enzyme storage buffer except substituting 25 mM imidazole for 25 mM EDTA). Protein purity was assessed using 4-12% gradient SDS-PAGE (Life Technologies) and total protein yield determined by A280.

### Enzyme assays and data analysis

For thermal stability assays, triplicate aliquots of freshly-purified proteins (less than 48 hours after cell lysis) in enzyme storage buffer at concentrations ranging from 0.01 to 0.1 mg/mL were incubated at 8 temperatures in 2.5 C increments between 30 C and 50 C in a thermal cycler (BioRad). After 30 minutes, proteins were transferred to assay plates containing 100 mM pNPG in enzyme storage buffer. The rate of pNP production was monitored over 60 minutes and product formation rate was determined by linear fit (Gen5).

For each BglB variant, triplicate rates at 8 temperatures were normalized to the [0,1] interval and fit to the logistic equation `1/(1+exp(-k*(x-x0)))` using least squares (accessed via SciPy's `curve_fit`) to determine the protein melting temperature, and 1 standard deviation error from fitted value. Michaelis-Menten parameters for each mutant were determined as described previously [Carlin]. Mutants for whom Michaelis-Menten parameters were previously reported use values from the referenced publication (assays were not repeated for all mutants).

# Results and discussion

## Soluble protein expression in the heterologous host E. coli

A total of 125 individually-designed mutations to the BglB sequence were produced, purified, and characterized in this study, along with 10 batched replicates of the native sequence. Of 125 mutant proteins synthesized and assayed, 92 express and purify as soluble protein. The remaining 33 mutants did not appear on a gel after at least 2 independent production attempts.

## Protein melting temperature for 69 mutants of BglB

Of the 92 solubly-expressed mutants, melting temperature for 69 was determined by fitting rate data collected at 8 temperatures from 30 to 50 C to the logistic equation. For the the remaining 21 mutants that expressed as soluble protein,  could not be used to determine Tm because the kinetic constants for these mutants are below our limit of detection.

The wild type BglB replicates had an average melting temperature of 39.6 C across assay data from 10 biological replicates. The average melting temperature was 39.0 ± 1.7 C, and all the Tms determined ranged from 32.9 to 45.3 C (12.4 C). Of 69 mutants for which Tm was determined, Tm for 37 of them falls within 1 degree C of the wild type Tm (58%). Of the remaining 49 Tm values, 18 exhibited a lower melting temperature than that of native BglB and 7 displayed a higher melting temperature. The highest Tm observed in this data set is N404C, which increased the Tm of BglB to 45.3 C (dTm=+5.4 ˚C), while the lowest Tm observed was for mutant H178A, which had a delta Tm of -7.0 C. **Figure 2** depicts the data set as a heat map.

A novel finding was a nearly three-degree increase in melting temperature by single point mutant N404C. The BglB crystal structure reveals a weak hydrogen bond between N404 and the backbone at L402. Molecular modeling of N404C predicts the loss of this hydrogen bond to the protein's alpha helix, in which the protein is allowed to repack into more energetically favorable states.

Point mutation E222H had a melting temperature of 34.7 C, a nearly five degree decrease than that of native BglB. Previous studies show strong hydrogen bond interaction, 2.6 and 3.1 Å, between E222 and its neighboring R240 residue [cite]. The introduction of histidine at the mutation site causes the loss of these strong hydrogen bonds as well as the creation of electrostatic repulsion between the partially positively charged and positively charged amino acids. The cumulative effect of this mutation results in the protein's decreased stability at lower temperatures.

Similarly, the point mutation W120F resulted in a delta Tm of +1.6 C. The BglB crystal structure indicates a weak hydrogen bond between W120 and the backbone of N163, which could be construed as an unsatisfied polar interaction. The mutation to Phe maintains the structural integrity at the mutation site as well as removes the unsatisfied hydrogen bond donor to the neighboring alpha helix. The increased stability is then due to destabilization of the unfolded state, which exposes a hydrophobic Phe to bulk solvent. It is also worth noting that the mutant W120A results in no soluble protein production after 3 independent attempts, indicating that W120 plays a key role in stabilizing the protein. Previous studies have shown a similar increase in stability upon W to Y point mutations [Fulton]. Analysis of a multiple sequence alignment of 1,554 proteins from the Pfam database reveals approximately equal probability of finding Trp, Phe, or Tyr here, and less than 1% representation of all of the other amino acids combined. Thus, the experimental measurements and the evolutionary record agree that W120 plays a key role in stabilizing BglB. No major structural changes are predicted via Rosetta modeling.  

## Kinetic constants for additional 25 variants of BglB expand previous data set

In addition to expression and thermal stability data, we also report kinetic constants for 12 mutants whose kinetic constants had not previously been determined. Together, data for soluble expression, kinetic constants (kcat, km, kcat/km), and melting temperature are reported for 128 mutants of BglB (Supporting Information). **Figure 2** depicts the data set as a heat map, with values relative to native BglB. Based on the maximum concentration of enzyme used in our assays and colorimetric absorbance changes at the highest substrate concentration used, we estimate our limit of detection for kcat/KM to be 10 M<sup>-1</sup>min<sup>-1</sup>.

![](figure_2.png)

**Figure 2: heat map showing the expression, Tm, kcat, KM, kcat/KM, and percent conservation for 128 mutants**. For expression (first column), a black box indicates soluble expression, and a white box indicates no expression for this mutant. The last column indicates the degree of conservation for the BglB native residue based on an alignment of 1,554 sequences from Pfam. For all other values, the heat map depicts each value relative to the wild type value. For Tm, a linear scale from -6 (purple) to 6 (yellow) is used to depict deltaTm compared to wild type. For kcat, KM, and kcat/km, a diverging colormap is used, with blue values indicating lower values, and yellow indicating higher values, as indicated by the color legend. For KM, the value 1/KM is plotted, so that a higher plotted value indicates greater affinity between BglB and pNPG for the mutation.

![](figure_3.png)

**Figure 3: Single point mutants with extreme effects on thermal stability**:  Four mutant panels are shown, sorted from left to right by delta-Tm. For each mutant, PyMOL images depict the local area of the mutation in the BglB WT protein (top) and Rosetta model of mutatation (bottom). In the bottom most panel, the assay for the mutant (black dots) is shown along with the fitted Tm (black dashed line) and the melting curve for WT BglB (green) as measured over 8 replicates in this study (the Y axis is in arbitrary units and the X axis covers the range 30-50 C).

## Data on a variety of mutants with preturbed eletric fields highlights the electrostatic weaknesses of current modeling protocols

Additionally, we have additional information about R240. The mutation R240E, which plates a glutamate near the existing glumatate at E222 that forms two hydrogen bonds to R240, placing two negatively-charged residues in close proximity, results in no soluble protein under our experimental conditions after at least 2 repeated attempts to produce and purify the protein, providing evidence for the belief that the mutation R240E is destabilizing. Whereas Rosetta modeling predicts no change in deltaG (within 1 REU) for this mutation, and example of a mutation that would have been missed using Rosetta. However, FoldX assigns a ddG of 3.78, highlighting the importance of a diverse set of estimates chosen by machine learning on experimental data to train algorithms. (Mutation R240D, which places a carboxylate only 1.54 A further away, results in a soluble protein.)

## Relationship between Tm and kcat, KM, and kcat/KM in BglB

In the BglB system, Tm and each of kcat, KM, and kcat/KM were not found to be correlated. Statistical analyses reveal that correlation between Tm and kcat (PCC -0.04), KM (-0.14), and kcat/KM (-0.02).

Similarly, in the case of conservation within the Pfam family 1 glycoside hydrolases, the correlation between percent conservation, as assessed by PCC, was found to be -0.275715, 0.202354, -0.279241, for kcat, KM, and kcat/KM, respecticly. For thermal stability and E. coli expression, the PCC values were 0.262210 and -0.239749.

![](figure_4.png)

**Figure 4: statistical independence of protein melting temperature (Tm) and kinetic constants.** Plots of Tm versus kinetic constant for each of kcat, KM, and kcat/km. Tm values are on a linear scale (C) and kinetic constant values are on a log scale.

## Function-stability relationships in BglB uncovered by mutations

The concept of stability-function tradeoff is a debated idea that some think to be a general governing principle of functional proteins such as enzymes. This idea stems from the observation that enzyme active sites, such as the one in BglB, often position charged residues in close proximity in order to present optimal geometry for catalysis. [Fersht] However, the putative function-stability tradeoffs of residues not directly involved in the reaction chemistry has been explored only at a small number (n<5) of sequence positions in T4 lysozyme [Matthews] and barnase [Fersht].

The data set reported here suggests a more nuanced view of the function-stability tradeoffs in BglB. Besides providing a way to assess the, statistical correlation between protein stability and kinetic constants, the individual mutants reported here reveal unexpected structure-function-stability relationships in BglB.

The mutation N404C reveals a clear tradeoff between function and stability: the mutation increases the functional Tm of BglB by 2.7 C while decreasing kcat by 10-fold. Similarly, the mutation W325L increases functional Tm by 2.6 C, while decreasing kcat by 9-fold. However, the opposite is also true in the BglB system. The mutation R240A does not change Tm, by increases kcat by 10-fold. For the Michaelis constant KM, a decrease (indicating greater enzyme-substrate affinity) of 10-fold in mutant N220Y is accompanied by an increase in Tm by 1.9 C.

Furthermore, the percentage of mutants in our data set that exhibit function-stability tradeoffs is 18.8% for kcat, 25.8% for KM, and 21.0% for kcat/KM, lower than the 50% that might be expected from chance alone. Since each mutation can be ++, --, or +-, -+ (both indicating a tradeoff), you would expect 50% tradeoffs if the effect of each mutation is independent.  

These differences may reflfect biases in the underlying data sets. For example, in the t4 data set, only about 18 single mutants are reported, and only 3 are desatbilizing. In contrast, in our data set we have about 100X more data, and the mutations are more balanced. Furthermore, in the T4 data set, the only stabilizing mutations E11M and E11A were found via X-ray crystallopgraphy to induce large conformational changes in the T4 lysozyme structure, making their stabiliting effect versus their functional effect difficult to interpret.  

# Conclusion

The concept of stability-function tradeoff is a widely-accepted idea that is thought to be a general governing principle of functional proteins such as enzymes. The data set that forms the basis of this study is the largest data set of enzyme single amino acid substitutions for which the kinetic constants kcat, KM, and kcat/KM, as well as protein thermal stability and soluble protein expression in *E. coli* have been explicitly measured for each mutation, enabling direct statistical assessment of the relationship between function and stability in BglB. In addition, the data set presented here uncovered unexpected structure-stability relationships in BglB, including mutants with significant effects on thermal stability not predicted by computational modeling and having no evolutionary precedent in available databases.

---

### Supporting information

  1. Table of protein expression, melting temperature, kinetic constants, and statistical analysis for each of 125 mutants (CSV)
  1. Michaelis-Menten plots for each mutant for which kinetic constants are reported (ZIP file containing images)
  1. Plot of protein melting for each mutant for which Tm is reported (ZIP file containing images)
  1. Images of gel band for each protein used in this study (ZIP file or single TIFF image)
  1. Additional information about experimental procedures (text)

---

## References

[antibody_corr] High-throughput measurement, correlation analysis, and machine-learning predictions for pH and thermal stabilities of Pfizer-generated antibodies

[Masso] Accurate prediction of stability changes in protein mutants by combining machine learning with structure based computational mutagenesis

[Berliner] http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4170975/

[PPI dataset] http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1003592

[molec_mech] Molecular mechanisms of disease-causing mutations JBC

[Beadle] Structural Bases of Stability–function Tradeoffs in Enzymes

[Tokuriki] How Protein Stability and New Functions Trade Off

[RNase paper]


## Unused stuff

We find that, perhaps as a consequence of the lack of large data sets quantitatively relating enzyme mutations to function, the prediction of valeus from our data set using state-of-the-art molecular modeling achieves a Pearson's R value of -0.49 (FoldX PSSM). FoldX is molecular modeling software that sues a scoring fuction that has been trained on and optimized for point mutations. Emperical force fields for point mutations to evaluate the correlation between computational predictions of protein stability and the measured stability of 128 mutants. Dozens of software packages are available for prediting the stability effects of mutations [molec_mech]. Some use features derived from evolutionary information and sequence alone [SIFT], while some incorporate structural information [PolyPhen, AUTO-MUTE, PopMuSiC]. Some, like FoldX, use force-field–based molecular modeling combined with a statistical potential that has been trained on experimental data for point mutations [FoldX]. Other previous work has incorporated machine learning algorithms [Masso, Berliner], using the metrics generated by the various software packages as features. None, however, has been independently evaluated using an independent, curated training set. No study has yet combined both a standardized approach to experimental characterization, and produced enough mutants to allow a sufficiently large training set to evaluate predictive ability.

## Computational modeling and feature generation using Rosetta and FoldX

Three separate approaches to modeling the 128 mutants were taken in this study. In the first, 100 molecular models of each mutant enzyme were generated using the Rosetta Molecular Modeling Suite by Monte Carlo optimization of total system energy using enzyme design protocols as reported in our previous work. The lowest 10 of total system energy were selected for feature generation and a total of 50 features were calculated for each protein model. In a second approach, mutants were generated and scored using benchmarked [Kellogg] Rosetta ddG application, which generates 50 models and averages the best 3 models to provide a feature set of 13 features. In a third approach, models were generated and scored using FoldX PSSM application, which generates a feature set of 13 features per model. The lowest energy model was selected for feature generation.

### Machine learning methods

Elastic net for Tm

Random forests for expression

## Sequence-stability-function relationships in BglB

In agreement with previous studies, mutants that did not express followed rules broadly consistent with previous work and sequence conservation, such as the large destabilizing effect of any substitution at W407. Other mutants that did not express mostly followed well-established observations of destabilizing effects, such as the introduction of proline into an alpha helix (Y166P, Q19P), the mutation of topology-defining/helix-ending proline residues (P329N), mutations to and from glycine (), the introduction of charged residues into the hydrophobic core (A236E, F72H), and extreme small-to-large mutations (A227W).

It is interesting to note that all the mutations we made in our study to the catalytic nucleophile (E353) and acid-base (E164) polar residues resulted in soluble protein, consistent with the idea that enzymes must trade structural stability for function in placing these two negatively-charged groups less than 5 Å apart (the sole exception, E164G, could be the result of an altered folding trajectory due to the conformational freedom of glycine). Notable is that all the catalytic knockouts, E164A, E164R, E353A, Y295A, Y295G, except for E164G, where a glycine is inserted in an alpha helix, result in protein that is solubly expressed in our system, consistent with the hypothesis that mutants of catalytic residues in the BglB system are not destabilizing.

http://www.pnas.org/content/92/2/452.full.pdf

Together, the importance of enzymes in biotechnology and human disease make the accurate modeling of enzyme stability an important goal of the protein modeling community. Accurate prediction of a point mutation's effect on enzyme stability would unlock rational protein engineering approaches, where the information could be used immediately to rationally engineer an enzyme's functional envelope for a desired situation as has been previously explored \cite{22575958}. Furthermore, understanding the changes in enzyme stability that occur upon point mutations would provide huge insight into understanding inherited diseases of metabolism [cite], cancer [cite], as well as mechanisms by which bacteria become resistant to antibiotics, which is lately a public health menace [cite].

Traditionally, studies of this type have sought to design more thermostable mutants. However, we are interested more in less thermostable mutants in the context of human disease. Many inherited disorders of metabolism are caused by random point mutations in germline protein-coding genes. Disease results when the mutant protein lacks a key catalytic residue, or when protein misfolding or instability leads to an inactive form. Since two alleles of metabolic enzymes are simultaneously expressed, a precise quantitative measurement of the kinetic activity and protein stabilty is necessary for the predicition to have any relevance.

The finding that protein melting temperatures and kinetic constants reported in this study were not found to be correlated. This sheds doubt on studies [Romero] that convolute these parameters into specific activity measurements, for one can never be sure the functional parameter that is tuned to result in an inactive mutant. In biotechnology, however, this finding is a boon: it means that kinetics and thermal stability, and perhaps other parameters of the functional envelope, arise from independent (if overlapping) biophysical properties, and thus can be rationally modulated independently. This neatly avoids the "two-objective optimisation problem" assumed to exist by engineers seeing to maximize two parameters (such as kcat and thermal stability) independently [http://pubs.rsc.org/en/content/articlepdf/2015/cs/c4cs00351a]. Whether this finding is relevant to BglB alone is unclear, as other studies have shown a strict linear correlation between delta-G and kcat/km for RNase [http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1219248/]. However, it will be impossible to know without large data sets such as the one presented here.
