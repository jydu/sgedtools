
Ungroup sites:
==============

```bash
python ../src/sged-ungroup.py \
         --sged=HOG000003295_charge_stats_pvalues.csv \
         --data=Size,IsConstant \
         --output=HOG000003295_sites.tsv
```

Create a PDB index:
===================

```bash
python ../src/sged-create-structure-index.py \
         --pdb=3BXY.pdb \
         --alignment=HOG000003295_bppalnscore.mase \
         --output=HOG000003295_PdbIndex.txt
```

Translate according to PDB:
===========================

```bash
python ../src/sged-translate-coords.py \
         --sged=HOG000003295_charge_stats_pvalues.csv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_charge_stats_pvalues_PDB.csv \
         --name=PDB

python ../src/sged-translate-coords.py \
         --sged=HOG000003295_sites.tsv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_sites_PDB.tsv \
         --name=PDB
```

Structural statistics:
======================

Compute 3D distances, RSA, and secondary structures:

```bash
python ../src/sged-pdb-infos.py \
         --sged=HOG000003295_charge_stats_pvalues_PDB.csv \
         --group=PDB \
         --pdb=3BXY.pdb \
         --chain=A \
         --measures=AlphaDist,DSSPsum \
         --output=HOG000003295_charge_stats_pvalues_PDB_infos.csv

```

RSA and structure per site:
```bash
python ../src/sged-pdb-infos.py \
         --sged=HOG000003295_sites_PDB.tsv \
         --group=PDB \
         --pdb=3BXY.pdb \
         --chain=A \
         --measures=DSSP \
         --output=HOG000003295_sites_PDB_infos.csv

``
