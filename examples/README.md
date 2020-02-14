
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
python ../src/sged-groups2structure.py \
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
