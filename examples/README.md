
Create a PDB index:
===================

```bash
python ../src/sged-groups2structure.py --pdb=3BXY.pdb --alignment=HOG000003295_bppalnscore.mase --output=HOG000003295_PdbIndex.txt
```

Ungroup sites:
==============

```bash
python ../src/sged-ungroup.py --sged=HOG000003295_charge_stats_pvalues.csv  --data=Size,IsConstant --output=HOG000003295_sites.tsv
```


