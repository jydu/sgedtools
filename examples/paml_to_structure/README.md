# Mapping Positively Selected Sites on a Protein Structure


The identification of positively selected sites on a protein sequence gives insights into
adaptive changes that have occurred over time, thus the evolutionary history of species
and the function of the genes or proteins. After identifying the positively selected sites
on a protein sequence, the next step is to map these positions onto a three-dimensional
structure of the protein. This provides insight into how amino acid changes and the physicochemical
properties - such as size, charge, and polarity - of the protein interact with each other,
and how these changes could potentially impact the protein's function. In this annotated
example, we used lysozyme sequences to detect positively selected sites using Codeml from
[PAML][1] (for Phylogenetic Analysis by Maximum Likelihood) by Yang (1997) version 4.9j
(February 2020) and mapped them onto a lysozyme protein structure.

## Inferring positively selected sites in the lysozyme gene c sequences

To infer the positively selected sites in the lysozyme gene sequences, we duplicated
the results of the analysis by Yang and Nielsen (2002). Here we used the lysozyme data
set of [Messier and Stewart (1997)][2] and files for the lysozyme gene sequences of 19 distinct
primate species which is referred to as the ‘‘large’’ data set by [Yang and Nielsen (2002)][3]
to infer the positively selected sites under the branch-site model. 
 
As suggested in the example file in PAML, we used branch-site model A to construct 
branch-site test 2 (branch-site test of positive selection) and tested the alternative
hypothesis which was already specified in the default control file (`lysozymeLarge.ctl`).
The identified positively selected sites were given under the NEB (Naive Empirical Bayesian)
and BEB (Bayes Empirical Bayes) [(Yang et al. 2005)][4] procedure at the end of the analysis.

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

The converting process is done with the sged-paml2sged.py script:

```bash
python3 sged-paml2sged.py -p mlc.txt -o ps.csv -m bayesian
```
where

* -p refers to the PAML result file,
* -o is the output file, which will be written in SGED format (TSV or CSV),
* -m statistical method for positively selected sites, either "naive" or "bayesian".



## Getting a PDB file and creating a PDB index

To map the identified positively selected sites onto a protein structure, the protein's
three-dimensional(3D) structure is required. The `sged-create-structure-index.py` program can
be used to find the best-matched protein structure and create a PDB index.
The program can read two different file formats for protein structures, Protein Data Bank (PDB) file and macromolecular Crystallographic Information File (mmCIF).

We use one of the sequences of the foreground (the lineage ancestral the colobine monkey) branch, 
which is "colobus_Cgu&Can", to find the best-matched protein structure, or in other words, the closest
protein structure to our sequence from the data bank. 
The nucleic acid sequence should be translated into amino acid sequence. In this example, the sequence 
of the foreground branch was inserted to the MEGA7 (MEGA7: Molecular Evolutionary Genetics Analysis version 7.0 for bigger datasets) [(Kumar, Stecher, and Tamura 2016)][5] 
from a text file, translated into the protein sequence, and exported as fasta format.

```bash
python3 sged-create-structure-index.py \
          -p *.pdb \
          -f PDB \
          -a colobus_aa.fas \
          -g fasta \
          -o PDB_index.txt \
          -x
```

The command lines arguments are the following:

* -p the input PDB file(s). Several files can be listed use multiple -p arguments, or a glob pattern can be specified. Alternatively: PDB ids can be specified when the format is set to `remote:`.
* -f the structure file format, one of PDB or MMCIF. When the `remote:` prefix is added (e.g. `-f remote:MMCIF`), structure files will be downloaded from the server.
* -a the alignment file for which the index should be created,
* -g the alignment file format. Any file format can be used supported by Bio.SeqIO, check [here](https://biopython.org/docs/1.76/api/Bio.SeqIO.html) : Clustal (.clustal, .clustalw, .aln), FASTA (.fasta, .fas, .fa, .fsa, .mpfa), NEXUS (.nexus, .nxs), IntelliGenetics (.ig, .mase)...
* -o the output file where to write the PDB index,
* -x only keep the chain(s) with the lowest amount of incomplete data (e.g. amino acids with missing lateral chain).

To create the index, the program will align all sequences from the input aligment with all chains from the selected PDB entries, and keep the best matching pair.

## Coordinates translastion

Lastly, identified positively selected sites should be translated according to the index using the `sged-translate-coords.py` program.
In this step, the converted SGED file and the PDB index file are needed. 

```bash
python3 sged-translate-coords.py \
          -s ps.csv \
          -o translated.csv \
          -i PDB_index.txt \
          -n PDB \
          -c
```

Command line arguments:

* -s the SGED input file,
* -o the output file that will include identified the positively selected sites and corresponding PDB coordinates,
* -i the PDB index file,
* -n the column name for corresponding PDB coordinates,
* -c tells the program to use CSV format instead of TSV for SGED files

An example of the output for our example (first 5 rows):
```
Group,PDB,amino_acid,probability
[14],[ARG14],R,0.859
[21],[ARG21],R,0.858
[23],[ILE23],I,0.853
[37],[GLY37],G,0.510
...
```


## References

[1]: Yang Z. PAML: a program package for phylogenetic analysis by maximum likelihood. 
Comput Appl Biosci. 1997 Oct;13(5):555-6. doi: 10.1093/bioinformatics/13.5.555. 
PMID: 9367129.

[2]: Messier, W., Stewart, CB. Episodic adaptive evolution of primate lysozymes. 
Nature 385, 151–154 (1997). https://doi.org/10.1038/385151a0

[3]: Ziheng Yang, Rasmus Nielsen, Codon-Substitution Models for Detecting Molecular Adaptation
at Individual Sites Along Specific Lineages, Molecular Biology and Evolution, Volume 19, 
Issue 6, June 2002, Pages 908–917

[4]: Ziheng Yang, Wendy S.W. Wong, Rasmus Nielsen, Bayes Empirical Bayes Inference of Amino Acid 
Sites Under Positive Selection, Molecular Biology and Evolution, Volume 22, Issue 4, April 2005, 
Pages 1107–1118

[5]: Kumar S, Stecher G, Tamura K. MEGA7: Molecular Evolutionary Genetics Analysis Version 7.0 for Bigger Datasets.
 Mol Biol Evol. 2016 Jul;33(7):1870-4. doi: 10.1093/molbev/msw054. Epub 2016 Mar 22. PMID: 27004904; PMCID: PMC8210823.
