# Investigation into the computational prediction of catalysis in a diverse set of oxidoreductase enzymes

## Introduction

While predicting the fold of a protein from sequence alone has proven an elusive goal for biochemistry, biophysics, and computer scientists working on molecular modeling, and even more challenging question is the prediction of functional enzyme-substrate pairs. As protein folding algorithms, led by the Rosetta Molecular Modeling Suite, grow more able to predict the folds of novel sequences [13 new folds by adding evolutionary information], we asked whether it is possible to assess simple, physically-realistic calculations based on Rosetta models of enzymes and correlate these to function on arbitrary substrates.

In order to begin to address this question, we looked at oxidoreductases, a family of enzymes that use NADH or NADPH as a cofactor in the reduction of double bonds in small molecules in the classes of alkenes, ketones, aldehydes, and imines. With the similar chemistry of hydride donation as a backdrop, we investigate the set of mechanisms used by this family of enzymes on a model panel of small substrates that were predicted to be compatible with the active sites of 30 oxidoreductase enzymes that form a diverse sample from the tree of life.

**Figure 1. A) Tree of life showing the enzymes selected for this study as labeled red nodes.** The tree of life was drawn from a recent Nature Microbiology paper. B) Overlay of the variety of NADH binding conformations in this set of enzymes. C) Overview of the structural diversity within the chosen sample set.

## Methods

### Diverse and representative set selection

In order to characterize a diverse set of proteins, we began with all structures in the PDB containing a bound NAD cofactor, that had been reported to be expressed in *E. coli*, and that were between 200 and 600 residues in length. Additionally, we filtered out P450-like and DHFR-like proteins since they are massively overrepresented in existing data sets. Of the 200-some structures meeting these critera in the PDB in March 2014, we selected 40 that had been reported previously to catalyze a hydride transfer from NAD(P)H to a small molecule substrate, based on reading each crystal structure paper and analyzing the mechanism.

### Molecular cloning

For each of 30 genes that were selected based on the above criteria, a gene was codon-optimized for *Escherichia coli* and synthesized as linear dsDNA (Life Technologies). The constructs were individually cloned into a pET29b+ vector linearized with XhoI and NdeI (NEB) using Gibson assembly master mix (NEB) and verified by Sanger sequencing (Operon).

### Protein production and purification

Sequence-verified constructs were transformed into chemically competent *Escherichia coli* BLR. Single colonies were used to inoculate 25 mL Terrific Broth (TB) (Fisher) in 50 mL Falcon tubes (Fisher) with breathable tube seals (Fisher). After incubation at 37 C with shaking at 300 RPM for 24 hours, cells were pelleted and media replaced with 25 mL autoinduction medium (recipe). After incubation at 18 C with shaking at 300 RPM for 24 hours, cells were pelleted, resuspended in enzyme storage buffer (50 mM HEPES, 150 mM sodium chloride, 25 mM EDTA, pH 7.50) and lysed with BugBuster protein extraction reagent (Novogen).

His-tagged proteins were individually purified via immobilized metal ion affinity chromatography using 50 µL bed volume of HisPur cobalt resin (Thermo) and eluted in 300 µL enzyme storage buffer. Protein expression was assessed using 4-12% gradient SDS-PAGE (Life Technologies).

### Substrate panel selection

A substrate panel covering the functional groups alkene, aldehyde, ketone, and imine was selected by modeling each of the proteins in the set and obtaining a consensus active site pocket. This small pocket was used as a guide in choosing a diverse panel of substrates containing double bonds that were predicted to be stericly compatible with the consensus binding pocket. The substrates chosen were

### Enzyme activity assays

Fresh stocks of NAD(P)H (as appropriate for protein) were prepared before each assay. 25 µL cofactor solution, 25 µL substrate stock, and 50 µL purified enzyme were mixed in wells of a 96-well plate. After mixing, the absorbance at X was monitored at 1 minute intervals for 60 minutes, and a line of best-fit was determined by BioTex Gen5 software

### DFT hydride affinity calculations

### Molecular modeling with Rosetta

## Results

### Protein expression of scaffold set

### Enzyme activity across a panel of substrates
