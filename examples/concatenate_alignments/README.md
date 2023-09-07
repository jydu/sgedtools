Objective:

Concatenate several input alignment files and create a SGED file keeping track of the coordinates.

We can use the `sged-concatenate-alignments` program to do all operations.

```bash
python3 ../../src/sged-concatenate-alignments.py \
    --aln1 ../data/cox1.aln.fasta \
    --aln2 ../data/cox2.aln.fasta \
    --output-aln cox12.aln.fasta \
    --output-sged cox12.index.sged
```

This creates a concatenated fasta alignment (`cox12.aln.fasta`) and a file containing all relative and absolute positions:

```
Group	AlnPos	AlnIndex	AlnId
[1]	[1]	1	cox1.aln.fasta
[2]	[2]	1	cox1.aln.fasta
[3]	[3]	1	cox1.aln.fasta
[4]	[4]	1	cox1.aln.fasta
[5]	[5]	1	cox1.aln.fasta
...
[518]	[518]	1	cox1.aln.fasta
[519]	[519]	1	cox1.aln.fasta
[520]	[520]	1	cox1.aln.fasta
[521]	[521]	1	cox1.aln.fasta
[522]	[522]	1	cox1.aln.fasta
[523]	[1]	2	cox2.aln.fasta
[524]	[2]	2	cox2.aln.fasta
[525]	[3]	2	cox2.aln.fasta
[526]	[4]	2	cox2.aln.fasta
[527]	[5]	2	cox2.aln.fasta
...
[748]	[226]	2	cox2.aln.fasta
[749]	[227]	2	cox2.aln.fasta
[750]	[228]	2	cox2.aln.fasta
[751]	[229]	2	cox2.aln.fasta
[752]	[230]	2	cox2.aln.fasta
```

To concatenate more than 2 alignments, we need to provide a list of files:

```bash
echo "../data/cox1.aln.fasta" > aln.lst
echo "../data/cox2.aln.fasta" >> aln.lst
echo "../data/cox3.aln.fasta" >> aln.lst
```

Then run the script:

```bash
python3 ../../src/sged-concatenate-alignments.py \
    --aln-list aln.lst \
    --output-aln coxall.aln.fasta \
    --output-sged coxall.index.sged
```

The resulting index looks like:

```bash
Group	AlnPos	AlnIndex	AlnId
[1]	[1]	1	cox1.aln.fasta
[2]	[2]	1	cox1.aln.fasta
[3]	[3]	1	cox1.aln.fasta
...
[520]	[520]	1	cox1.aln.fasta
[521]	[521]	1	cox1.aln.fasta
[522]	[522]	1	cox1.aln.fasta
[523]	[1]	2	cox2.aln.fasta
[524]	[2]	2	cox2.aln.fasta
[525]	[3]	2	cox2.aln.fasta
...
[750]	[228]	2	cox2.aln.fasta
[751]	[229]	2	cox2.aln.fasta
[752]	[230]	2	cox2.aln.fasta
[753]	[1]	3	cox3.aln.fasta
[754]	[2]	3	cox3.aln.fasta
[755]	[3]	3	cox3.aln.fasta
...
[1011]	[259]	3	cox3.aln.fasta
[1012]	[260]	3	cox3.aln.fasta
[1013]	[261]	3	cox3.aln.fasta
```

Create a phylogenetic tree from the concatenated alignment. We use the PhyML program withing SeaView (LG+Gamma model, best of NNI + SPR, aLRT bootstraps), and save the unrooted tree. The alignment is pretty good so all sites are used. Using bppPhyView, we collapse all nodes with a bootstrap value lower then 0.6.
We then use comap to re-fit the model and perform substitution mapping in order to detect coevolving sites.

comap param=comap.bpp
Rscript computePValues.R

We create an index for each input alignments. We start with cox1, for which we download the structure 6J8M:

```bash
python3 ../../src/sged-create-structure-index.py \
         --pdb-id 6J8M \
         --pdb-format remote:mmCif \
         --alignment cox1.aln.fasta \
         --output cox1_PdbIndex.txt
```
This maps to chain A.
We then do the two other units, using the downloaded structure.

```bash
python3 ../../src/sged-create-structure-index.py \
         --pdb 6j8m.cif \
         --pdb-format mmCif \
         --alignment ../data/cox2.aln.fasta \
         --output cox2_PdbIndex.txt

python3 ../../src/sged-create-structure-index.py \
         --pdb 6j8m.cif \
         --pdb-format mmCif \
         --alignment ../data/cox3.aln.fasta \
         --output cox3_PdbIndex.txt
```

cox2 maps to chain B and cox3 to chain C.

We then combine the indexes to map each PDB chain to the concatenated alignment:

```bash
python3 ../../src/sged-liftover-index.py \
         --index1 cox1.aln.fasta_AlnIndex.txt \
         --index2 cox1_PdbIndex.txt \
         --output cox1.aln.fasta_PdbIndex.txt
```
TODO

We then merge the three indexes into a single one, since they are complementary:

TODO

Finally, we get the PDB coordinates of each coevolving group:

TODO

Are there inter subunit coevolution? show them on the structure...
