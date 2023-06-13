# RASER to structure

The RAte Shift EstimatoR [RASER](https://www.tau.ac.il/~penn/raser.html)[^1] is a likelihood-based method to test and infer site-specific evolutionary rate shifts. 
In some cases, it can be interesting to combine this analysis with a three dimensional structure analysis of proteins to detect site-specific changes in specific protein structures. 
Therefore, you may want to convert your RASER outputs into SGED files.

To do that, there are several steps:

## Convert Raser output to an SGED file

The first step is to transform the raser output file into a SGED package -readable file. It is done using the `sged-raser2sged.py` script.

```{bash}
python3 sged-raser2sged.py \
        -r CLU_000422_3_3_results_file.txt \
        -a CLU_000422_3_3.mase \
        -f ig \
        -o CLU_000422_3_3_sged.csv \
        -c
```
The command line arguments are:
* '-r' or '--raser=' is the raser result file.
* '-a' or '--alignment=' is the multiple seqeunce alignment file used to run RASER.
* '-f' or '--alignment-fomat=' is the format of the multiple sequence alignment file. It can be any file format supported by BioSeqIO ([check here](https://biopython.org/docs/1.76/api/Bio.SeqIO.html)): CLUSTAL (.clustal, .clustalw, .aln), FASTA (.fasta, .fas, .fa, .fsa, .mpfa), NEXUS (.nexus, .nxs), IntelliGenetics (.ig, .mase), ...
* '-o' or '--output=' is the name of the output file (SGED file)
* '-c' or '--csv'. This argument is optional. It's use will print the output file in CSV format. The default value is TSV.

This should make an output file looking like this (in TSV format): (several lines in the file)
```
Group   amino_acid      probability     Proba>0.95
[0]     -       NA      NA
[1]     -       NA      NA
...
[42]    G       0.54    
[43]    E       0.96    *
...
[590]   P       0.79    
[591]   D       0.78    
```

## Getting a PDB file and creating a PDB index

Once the SGED file is made, as we want to map the site-specific rate shifts onto the protein structure, we have to get the three dimensional structure of the protein. This 3D structure is given in files that can be fetched by the `sged-create-structure-index.py` program. This program can be used to find the best-matched protein structure and create a PDB index. The index is the best alignment of the chains in the protein structure file and the multiple sequence alignment. It can take several protein structure references and select the best match using two different file formats: Protein Data Bank (PDB) or macromolecular Crystallographic information file (mmCif).

```{bash}
python3 sged-create-structure-index.py \
        -p 5jca.pdb \
        -p 5jfc.pdb \
        -f PDB \
        -a CLU_000422_3_3.mase \
        -g ig \
        -o CLU_000422_3_3_index.txt \ 
        -x
```
The command line arguments are:
* '-p' or '--pdb=' is the name of the PDB file. Several files can be listed using multiple '-p' arguments or a global pattern can be specified (`-p pdbfile1 -p pdbfile2 ...`). Alternatively, PDB ids can be specified when the format is set to `remote:` (see example below).
* '-f' or '--format=' is the structure file format. It can be one of 'pdb' or 'mmCif'. If the `remote:` prefix is added, the struture files will be downloaded from the server.
* '-a' or '--alignment=' is the multiple sequence alignment file used to run RASER on which the index should be created.
* '-g' or '--alignment-format=' is the format of the multiple sequence alignment file. It can be any file format supported by BioSeqIO ([check here](https://biopython.org/docs/1.76/api/Bio.SeqIO.html)): CLUSTAL (.clustal, .clustalw, .aln), FASTA (.fasta, .fas, .fa, .fsa, .mpfa), NEXUS (.nexus, .nxs), IntelliGenetics (.ig, .mase), ...
* '-o' or '--output=' is the name of the output file (index file)
* '-x' or '--exclude-incomplete' only keeps the chain(s) with the lowest amount of incomplete data (e.g. amino acids with a missing lateral chain).

This example will be fetching the PDB files corresponding to the PDB ids on a specified server. It will then chose the best-matching protein sequence according to the alignment.
```{bash}
python3 sged-create-structure-index.py \
        -p 5JCA \
        -p 5JFC \
        -f remote:PDB \
        -a CLU_000422_3_3.mase \
        -g ig \
        -o CLU_000422_3_3_index.txt \ 
        -x
```

This step should create an index file looking like this: (first 10 lines)
```
# SGED index file version 0.99
# SGED input alignment = ../Database_test/Archaea/CLU_000422_3_3/07-Consensus/CLU_000422_3_3.mase
# SGED input alignment sequence = PYRFU_1.PE1332
# SGED input PDB = ./pdb5jca.ent
# SGED input PDB chain = L
# SGED index start
AlnPos,PdbRes
27,NA
28,PRO2
29,ARG3
...
```

## Coordinates translation

Once the index is created, we should translate the alignemnt positions selected in the previously created SGED file into chain sites in the index using the `sged-translate-coords.py`program.

```{bash}
python3 sged-translate-coords.py \
        -s CLU_000422_3_3_sged.csv \
        -i CLU_000422_3_3_index.txt \
        -n PDB \
        -o CLU_000422_3_3_translated_coords.csv \
        -c
```
The command line arguments are:
* '-s' or '--sged=' is the SGED input file comming from the raser conversion.
* '-i' or '--index=' is the index file previously created
* '-n' or '--name=' is the column name that will be given to the PDB coordinates.
* '-o' or '--output' is the name of the output file (translated coordinates file).
* '-c' or '--csv'. This argument is optional. It's use will print the output file in CSV format. The default value is TSV.

This should make an output file looking like this (in TSV format): (few lines)
```
Group	PDB	amino_acid	probability	Proba>0.95
[0]	[NA]	-	NA	NA
[1]	[NA]	-	NA	NA
[2]	[NA]	-	NA	NA
[3]	[NA]	-	NA	NA
[4]	[NA]	-	NA	NA
[5]	[NA]	-	NA	NA
[6]	[NA]	-	NA	NA
[7]	[NA]	-	NA	NA
[8]	[NA]	-	NA	NA
...
[41]	[SER15]	A	0.81	
[42]	[VAL16]	G	0.54	
[43]	[GLY17]	E	0.96	*
[44]	[GLU18]	R	0.76	
...
```

## Getting the Secondary Structure Information

The last step is to get the secondary structure information for the translated coordinates to be able to map the site specific evolutionary rate shifts. To do this, we will be using the `sged-structure-infos.py` script to calculate several informations on the 3D structure of the protein.

```{bash}
python3 sged sged-structure-infos.py \
        -s CLU_000422_3_3_translated_coords.csv \
        -p 5jca.pdb \
        -f PDB \
        -g PDB \
        -a L \
        -m DSSP \
        -o CLU_000422_3_3_structinfos
```

The command line arguments are:
* '-s' or '--sged=' is the SGED input file comming from the translation of coordinates.
* '-p' or '--pdb=' is the name of the PDB file to take as reference for the calculations. Of course, in this step too, you can use the `remote:` option just like in the previous program. The PDB reference to use is the one given in the beginning of the index file (line 4).
* '-f' or '--pdb-format=' is the structure file format. It can be one of 'pdb' or 'mmCif'. If the `remote:` prefix is added, the struture files will be downloaded from the server.
* '-g' or '--groups=' is the name of the column in which to get the index references from the translated coordinates file.
* '-a' or "--chain=' is the chain to consider in the PDB reference file. In the same way as the PDB reference to use, it can be found in the beginning of the index file (line 5).
* '-m' or '--measures=' is the names of the measures to do on the PDB file. It can be a list of the following items: **AlphaDist**, **ContactSubgraphs**, **ContactMap**, **DSSPsum**, DSSP, **ResidueDept** and SecondaryStructureLabels. (the bold measures are the ones that can only be done on groups so usually not used when converting a file from raser to a structure index. The SecondaryStructureLabels is to be used only if you have a `mmCif` reference file.). To do a list of these measures, type separate them using a comma, but no space, as shown here: `AlphaDist,DSSP,SecondaryStructureLabels`
* '-o' or '--output=' is the name of the output file containing the information on the structures.
* '-c' or '--csv'. This argument is optional. It's use will print the output file in CSV format. The default value is TSV. 

The ouptut file should look like this: (few lines)
```
Group	PDB	amino_acid	probability	Proba>0.95	Rsa	SecondaryStructure
[0]	[NA]	-	NA	NA	NA	NA
[1]	[NA]	-	NA	NA	NA	NA
[2]	[NA]	-	NA	NA	NA	NA
[3]	[NA]	-	NA	NA	NA	NA
[4]	[NA]	-	NA	NA	NA	NA
[5]	[NA]	-	NA	NA	NA	NA
[6]	[NA]	-	NA	NA	NA	NA
[7]	[NA]	-	NA	NA	NA	NA
[8]	[NA]	-	NA	NA	NA	NA
...
[29]	[ARG3]	N	0.6	NA	0.9475806451612904	-
[30]	[LEU4]	R	0.87	NA	0.22560975609756098	-
[31]	[ILE5]	V	0.83	NA	0.6390532544378699	-
[32]	[LYS6]	G	0.67	NA	0.48292682926829267	-
[33]	[ASP7]	A	0.83	NA	0.7607361963190185	S
[34]	[ARG8]	F	0.79	NA	0.3790322580645161	P
...
```

Once you have done this, you have the site-specific evolutionnary rate shift mapped on the 3D structure of a protein. 



## Citations

> [^1] Penn O, Stern A, Rubinstein ND, Dutheil J, Bacharach E, Galtier N., Pupko T. 2008. Evolutionary Modeling of Rate Shifts Reveals Specificity Determinants in HIV-1 Subtypes.
PLoS Computational Biology 4(11): e1000214
