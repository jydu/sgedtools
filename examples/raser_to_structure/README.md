# RASER to structure

The RAte Shift EstimatoR [RASER](https://www.tau.ac.il/~penn/raser.html)[^1] is a likelihood-based method to test and infer site-specific evolutionary rate shifts. 
In some cases, it can be interesting to combine this analysis with a three dimensional structure analysis of proteins to detect site-specific changes in specific protein structures. 
Therefore, you may want to convert your RASER outputs into SGED files.

To do that, there are several steps, which we illustrate with a sequence alignment of bacterial oxidoreductases. 

## Convert Raser output to an SGED file

The first step is to transform the raser output file into a SGED file. It is done using the `sged-raser2sged.py` script.

```{bash}
python3 ../../src/sged-raser2sged.py \
        -r CLU_000422_3_3_results_file.txt \
        -a CLU_000422_3_3.mase \
        -f ig \
        -o CLU_000422_3_3.sged \
```
('ig' is the python name for the Mase format.)
The program takes two input files:

* the raser result file.
* the multiple sequence alignment file used to run RASER.

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
python3 ../../src/sged-create-structure-index.py \
        -i 5JCA \
        -i 5JFC \
        -f remote:PDB \
        -a CLU_000422_3_3.mase \
        -g ig \
        -o CLU_000422_3_3_PdbIndex.txt \ 
        -x
```

This step should create an index file looking like this: (first 10 lines)

```
# SGED index file version 1.00
# SGED input alignment = CLU_000422_3_3.mase
# SGED input alignment sequence = PYRFU_1.PE1332
# SGED input PDB = ./pdb5jca.ent
# SGED input PDB chain = L
# SGED index start
AlnPos,PdbRes
27,NA
28,L:PRO2
29,L:ARG3
...
```

## Coordinates translation

Once the index is created, we translate the alignment positions selected in the previously created SGED file into chain sites in the index using the `sged-translate-coords.py`program.

```{bash}
python3 ../../src/sged-translate-coords.py \
        -s CLU_000422_3_3.sged \
        -i CLU_000422_3_3_PdbIndex.txt \
        -n PDB \
        -o CLU_000422_3_3_PDB.sged
```

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
python3 ../../src/sged-structure-infos.py \
        -s CLU_000422_3_3_PDB.sged \
        -p pdb5jca.ent \
        -f PDB \
        -g PDB \
        -m DSSP \
        -o CLU_000422_3_3_structinfos.sged
```

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
[29]	[ARG3]	N	0.6 NA	0.9475806451612904	-
[30]	[LEU4]	R	0.87	NA	0.22560975609756098	-
[31]	[ILE5]	V	0.83	NA	0.6390532544378699	-
[32]	[LYS6]	G	0.67	NA	0.48292682926829267	-
[33]	[ASP7]	A	0.83	NA	0.7607361963190185	S
[34]	[ARG8]	F	0.79	NA	0.3790322580645161	P
...
```

The site-specific evolutionnary rate shift are now mapped on the 3D structure of a protein. 
We can test whether secondary structure affects the probability of rate shifts. We compare sites with a probability of rate shift above 0.99:

```r
dat <- read.table("CLU_000422_3_3_structinfos.sged", header = T)
tbl <- table(dat$probability > 0.99, dat$SecondaryStructure)
chisq.test(tbl, simulate.p.value = TRUE)
```

```
> tbl
       
          -   B   E   G   H   I   P   S   T
  FALSE  81   7  72  13 149   3  11  37  42
  TRUE    1   2   4   0  12   2   0   3   4
```

```
	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  tbl
X-squared = 19.67, df = NA, p-value = 0.02049
```

There is an excess of rate shifts in secondare structure motifs in general.

## Citations

> [^1] Penn O, Stern A, Rubinstein ND, Dutheil J, Bacharach E, Galtier N., Pupko T. 2008. Evolutionary Modeling of Rate Shifts Reveals Specificity Determinants in HIV-1 Subtypes.
PLoS Computational Biology 4(11): e1000214
