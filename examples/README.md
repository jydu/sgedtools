
Ungroup sites:
==============

```bash
python ../src/sged-ungroup.py \
         --sged=HOG000003295_charge_stats_pvalues.csv \
         --data=Size \
         --output=HOG000003295_sites.tsv
```

List all residues in a PDB file:
================================

```bash
python ../src/sged-structure-list.py \
         --pdb=3EG4.pdb \
         --output=3EG4_sites.tsv \
         --chain=A
```

Create a PDB index:
===================

We compare all chains in the three available structures to find the best match in the alignment:
```bash
python3 ../src/sged-create-structure-index.py \
         --pdb=3BXY.pdb \
         --pdb=3EG4.pdb \
         --pdb=3GOS.pdb \
         --alignment=HOG000003295_bppalnscore.mase \
         --output=HOG000003295_PdbIndex.txt
```
which does the same as
```bash
python3 ../src/sged-create-structure-index.py \
         --pdb=*.pdb \
         --alignment=HOG000003295_bppalnscore.mase \
         --output=HOG000003295_PdbIndex.txt
```

Translate according to PDB:
===========================

```bash
python3 ../src/sged-translate-coords.py \
         --sged=HOG000003295_charge_stats_pvalues.csv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_charge_stats_pvalues_PDB.csv \
         --name=PDB

python3 ../src/sged-translate-coords.py \
         --sged=HOG000003295_sites.tsv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_sites_PDB.tsv \
         --name=PDB
```

Structural statistics:
======================

Compute 3D distances, RSA, and secondary structures:

```bash
python ../src/sged-structure-infos.py \
         --sged=HOG000003295_charge_stats_pvalues_PDB.csv \
         --group=PDB \
         --pdb=3BXY.pdb \
         --chain=A \
         --measures=AlphaDist,DSSPsum \
         --output=HOG000003295_charge_stats_pvalues_PDB_infos.tsv

```

RSA and structure per site:
```bash
python ../src/sged-translate-coords.py \
         --sged=HOG000003295_charge_siteInfos.csv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_charge_siteInfos_PDB.tsv \
         --name=PDB
python ../src/sged-structure-infos.py \
         --sged=HOG000003295_charge_siteInfos_PDB.tsv \
         --group=PDB \
         --pdb=3BXY.pdb \
         --chain=A \
         --measures=DSSP \
         --output=HOG000003295_charge_siteInfos_PDB_infos.tsv

``

Merging files:
==============

Get significant groups only:
```bash
csvgrep -t -c FDR -m yes HOG000003295_charge_stats_pvalues.csv | csvformat -T > HOG000003295_charge_significant.tsv
```

Get sites:
```bash
python ../src/sged-ungroup.py \
         --sged=HOG000003295_charge_significant.tsv \
         --data=Size \
         --output=HOG000003295_charge_significant_sites.tsv
```

```bash
python ../src/sged-merge.py \
         --sged1=HOG000003295_charge_siteInfos_PDB_infos.tsv \
         --sged2=HOG000003295_charge_significant_sites.tsv \
         --group=Group \
         --output=HOG000003295_charge_merged.tsv

```

Get distance between all pairs of sites in a PDB:
=================================================

```bash
python3 ../src/sged-structure-list.py \
        --pdb=3EG4.pdb \
        --chain=A \
        --output=3EG4_sites.tsv

python3 ../src/sged-get-all-pairs.py \
        --sged=3EG4_sites.tsv \
        --group=PDB \
        --output=3EG4_sitepairs.tsv

python3 ../src/sged-structure-infos.py \
        --sged=3EG4_sitepairs.tsv \
        --group=PDB \
        --pdb=3EG4.pdb \
        --chain=A \
        --measures=AlphaDist,DSSPsum \
        --output=3EG4_sitepairs_infos.tsv

```

Test if groups in one file are contained in groups from another:
================================================================

First add PDB coordinates to the coevolving groups:
```bash
python3 ../src/sged-translate-coords.py \
         --sged=HOG000003295_charge_significant.tsv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_charge_significant_PDB.tsv \
         --name=PDB
```

Then test if for all pairs in the 3D structure are in coevolving groups:
```bash
python3 ../src/sged-group-test-inclusion.py \
        --sged1=3EG4_sitepairs_infos.tsv \
        --group1=PDB \
        --sged2=HOG000003295_charge_significant_PDB.tsv \
        --group2=PDB \
        --output=3EG4_sitepairs_infos_coevolving.tsv \
        --result=IsCoevolving

```



Conversion from other software:
===============================

## DisEMBL

```bash
python3 ../src/sged-disembl2sged.py -d HOG000003295_scores.tsv -o HOG000003295_scores.sged
```
Convert to alignment positions. First create an index:

```bash
python3 ../src/sged-create-sequence-index.py -a HOG000003295_bppalnscore.mase -r seq451 -o HOG000003295_SeqIndex.txt
```
(seq451 was the sequence used to predict intrinsically disordered regions)

```bash
python3 ../src/sged-translate-coords.py \
         --sged=HOG000003295_scores.sged \
         --index=HOG000003295_SeqIndex.txt \
         --output=HOG000003295_scores_ref.sged \
         --name=AlnPos \
         --reverse
```


