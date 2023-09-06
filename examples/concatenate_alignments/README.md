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
[1]	[1]	1	../data/cox1.aln.fasta
[2]	[2]	1	../data/cox1.aln.fasta
[3]	[3]	1	../data/cox1.aln.fasta
[4]	[4]	1	../data/cox1.aln.fasta
[5]	[5]	1	../data/cox1.aln.fasta
...
[518]	[518]	1	../data/cox1.aln.fasta
[519]	[519]	1	../data/cox1.aln.fasta
[520]	[520]	1	../data/cox1.aln.fasta
[521]	[521]	1	../data/cox1.aln.fasta
[522]	[522]	1	../data/cox1.aln.fasta
[523]	[1]	2	../data/cox2.aln.fasta
[524]	[2]	2	../data/cox2.aln.fasta
[525]	[3]	2	../data/cox2.aln.fasta
[526]	[4]	2	../data/cox2.aln.fasta
[527]	[5]	2	../data/cox2.aln.fasta
...
[748]	[226]	2	../data/cox2.aln.fasta
[749]	[227]	2	../data/cox2.aln.fasta
[750]	[228]	2	../data/cox2.aln.fasta
[751]	[229]	2	../data/cox2.aln.fasta
[752]	[230]	2	../data/cox2.aln.fasta
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
[1]	[1]	1	../data/cox1.aln.fasta
[2]	[2]	1	../data/cox1.aln.fasta
[3]	[3]	1	../data/cox1.aln.fasta
...
[520]	[520]	1	../data/cox1.aln.fasta
[521]	[521]	1	../data/cox1.aln.fasta
[522]	[522]	1	../data/cox1.aln.fasta
[523]	[1]	2	../data/cox2.aln.fasta
[524]	[2]	2	../data/cox2.aln.fasta
[525]	[3]	2	../data/cox2.aln.fasta
...
[750]	[228]	2	../data/cox2.aln.fasta
[751]	[229]	2	../data/cox2.aln.fasta
[752]	[230]	2	../data/cox2.aln.fasta
[753]	[1]	3	../data/cox3.aln.fasta
[754]	[2]	3	../data/cox3.aln.fasta
[755]	[3]	3	../data/cox3.aln.fasta
...
[1011]	[259]	3	../data/cox3.aln.fasta
[1012]	[260]	3	../data/cox3.aln.fasta
[1013]	[261]	3	../data/cox3.aln.fasta
```


