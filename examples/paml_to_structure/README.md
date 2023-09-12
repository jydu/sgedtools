# Mapping Positively Selected Sites on a Protein Structure


The identification of positively selected sites on a protein sequence gives insights into
adaptive changes that have occurred over time, thus the evolutionary history of species
and the function of the genes or proteins. After identifying the positively selected sites
on a protein sequence, the next step is to map these positions onto a three-dimensional
structure of the protein. This provides insight into how amino acid changes and the physicochemical
properties - such as size, charge, and polarity - of the protein interact with each other,
and how these changes could potentially impact the protein's function. In this annotated
example, we used lysozyme sequences to detect positively selected sites using Codeml from
PAML[^1] (for Phylogenetic Analysis by Maximum Likelihood) by Yang (1997) version 4.9j
(February 2020) and mapped them onto a lysozyme protein structure.

## Inferring positively selected sites in the lysozyme gene c sequences

To infer the positively selected sites in the lysozyme gene sequences, we duplicated
the results of the analysis by Yang and Nielsen (2002)[^3]. Here we used the lysozyme data
set of Messier and Stewart (1997)[^2] and files for the lysozyme gene sequences of 19 distinct
primate species which is referred to as the ‘‘large’’ data set by Yang and Nielsen (2002)
to infer the positively selected sites under the branch-site model. 
 
As suggested in the example file in PAML, we used branch-site model A to construct 
branch-site test 2 (branch-site test of positive selection) and tested the alternative
hypothesis which was already specified in the default control file (`lysozymeLarge.ctl`).
The identified positively selected sites were given under the NEB (Naive Empirical Bayesian)
and BEB (Bayes Empirical Bayes) (Yang et al. 2005)[^4] procedure at the end of the analysis.

The ancestor of the colobine monkeys was used as the foreground branch and the rest in the phylogeny were
background branches as it was in the original analysis. Please check the tree file
"lysozymeLarge.trees" in the lysozyme example in PAML.

The PAML user's guide can be found [here](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf).

After installing PAML, the analysis can be run using the default control file for the large dataset:

```bash
codeml lysozymeLarge.ctl
```


## Converting PAML output to the SGED format

To map the positively selected sites onto a protein structure, the identified positively
selected sites in the main result file of PAML should be extracted and converted to the
SGED format.
In this example, the SGED file will contain all (positively selected) sites reported by PAML.

The converting process is done with the `sged-paml2sged.py` script:

```bash
python3 ../../src/sged-paml2sged.py \
    --paml mlc \
    --output lysozymeLarge-possel.sged \
    --method bayesian
```
We retrieve the sites detected by the empirical Bayesian method.


## Getting a PDB file and creating a PDB index

To map the identified positively selected sites onto a protein structure, the protein's
three-dimensional(3D) structure is required. The `sged-create-structure-index.py` program can
be used to find the best-matched protein structure and create a PDB index.
The program can read two different file formats for protein structures, Protein Data Bank (PDB) file and macromolecular Crystallographic Information File (mmCIF).

We use one of the sequences of the foreground (the lineage ancestral the colobine monkey) branch, 
which is "colobus_Cgu&Can", to find the best-matched protein structure, or in other words, the closest
protein structure to our sequence from the data bank. 
The nucleic acid sequence should be translated into amino acid sequence. In this example, the sequence 
of the foreground branch was inserted to the MEGA7 (MEGA7: Molecular Evolutionary Genetics Analysis version 7.0 for bigger datasets; Kumar, Stecher, and Tamura 2016)[^5] 
from a text file, translated into the protein sequence, and exported as fasta format.

```bash
python3 ../../src/sged-create-structure-index.py \
          --pdb "*.pdb" \
          --pdb-format PDB \
          --alignment colobus_aa.fas \
          --alignment-format fasta \
          --output lysozymeLarge_PdbIndex.txt \
          --exclude-incomplete
```

Note that we use a glob pattern to analyse all PDB files in the folder. For the glob pattern to be parsed by the program and not bash (which would lead to a program error), we surround the pattern by quotes.
To create the index, the program will align all sequences from the input aligment with all chains from the selected PDB entries, and keep the best matching pair.

## Coordinates translation

Lastly, identified positively selected sites should be translated according to the index using the `sged-translate-coords.py` program.
In this step, the converted SGED file and the PDB index file are needed. 

```bash
python3 ../../src/sged-translate-coords.py \
          --sged lysozymeLarge-possel.sged \
          --output lysozymeLarge-possel-PDB.sged \
          --index lysozymeLarge_PdbIndex.txt \
          --name PDB
```

Which gives:

```
roup   PDB     amino_acid      probability
[14]    [A:ARG14]       R       0.859
[21]    [A:ARG21]       R       0.858
[23]    [A:ILE23]       I       0.853
[37]    [A:GLY37]       G       0.510
[41]    [NA]    R       0.710
[50]    [NA]    R       0.704
[62]    [A:ARG62]       R       0.564
[87]    [A:ASP87]       D       0.869
[126]   [A:GLN126]      Q       0.710
```

## Visualizing the 3D Protein Structure

We can visualize the 3D protein structure and positively selected sites onto it
using PyMOL (Schrödinger et al. 2020)[^6] version 2.5.5 (2023).

First, we open the PDB file: 
```bash
pymol /path/to/134l.pdb
```

The model and color of the protein structure can be arranged by using options in PyMOL or using the PyMOL command line interface:

```
show mesh, all
```

Then we select the positively selected residues to show onto the protein structure
```
select my_residues, resi 14+21+23+37...
```

Again the model and color of the positively selected residues can be arranged by PyMOL options or,
```
show spheres, my_residues
```

Finally, the 3D structure can be saved in PNG format or PDB format
```
png /path/to/final.png, dpi=300, ray=1
save /path/to/final.pdb, selection=my_residues
```

```

[^1]: Yang Z. PAML: a program package for phylogenetic analysis by maximum likelihood. 
Comput Appl Biosci. 1997 Oct;13(5):555-6. doi: 10.1093/bioinformatics/13.5.555. 
PMID: 9367129.

[^2]: Messier W and Stewart CB. Episodic adaptive evolution of primate lysozymes. 
Nature 385, 151–154 (1997). https://doi.org/10.1038/385151a0

[^3]: Yang Z and Nielsen R. Codon-Substitution Models for Detecting Molecular Adaptation
at Individual Sites Along Specific Lineages, Molecular Biology and Evolution, Volume 19, 
Issue 6, June 2002, Pages 908–917

[^4]: Yang Z, Wong WSW, and Nielsen R. Bayes Empirical Bayes Inference of Amino Acid 
Sites Under Positive Selection, Molecular Biology and Evolution, Volume 22, Issue 4, April 2005, 
Pages 1107–1118

[^5]: Kumar S, Stecher G, and Tamura K. MEGA7: Molecular Evolutionary Genetics Analysis Version 7.0 for Bigger Datasets.
 Mol Biol Evol. 2016 Jul;33(7):1870-4. doi: 10.1093/molbev/msw054. Epub 2016 Mar 22. PMID: 27004904; PMCID: PMC8210823.
 
[^6]: Schrödinger L and DeLano, W. (2020). PyMOL. Retrieved from http://www.pymol.org/pymol
 
