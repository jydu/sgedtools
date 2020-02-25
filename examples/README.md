
Ungroup sites:
==============

```bash
python ../src/sged-ungroup.py \
         --sged=HOG000003295_charge_stats_pvalues.csv \
         --data=Size \
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
         --output=HOG000003295_charge_stats_pvalues_PDB_infos.tsv

```

RSA and structure per site:
```bash
python ../src/sged-translate-coords.py \
         --sged=HOG000003295_charge_siteInfos.csv \
         --index=HOG000003295_PdbIndex.txt \
         --output=HOG000003295_charge_siteInfos_PDB.tsv \
         --name=PDB
python ../src/sged-pdb-infos.py \
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
