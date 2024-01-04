# Analysis of the evolutionary rate of intrinsically disordered regions

This example demonstrates how to use the sgetools to analyse the output of the DISEMBL program to predict intrinsically disordered regions in a protein.
We use the HOG000003295 family from [^1].

We assume that all dependencies have been installed into a conda environment and activate it:

```bash
conda activate sgedtools-env
```

## Convert to SGED format

```bash
sged-disembl2sged \
    --disembl ../data/HOG000003295_scores.tsv \
    --output HOG000003295_scores.sged
```

## Convert to alignment positions

First create an index:

```bash
sged-create-sequence-index \
    --alignment ../data/HOG000003295_bppalnscore.mase \
    --alignment-format ig \
    --reference seq451 \
    --output HOG000003295_SeqIndex.txt
```
(seq451 was the sequence used to predict intrinsically disordered regions.)
Then use it to translate positions. We use the index in reverse mode, as we want to translate from a specific sequence to alignment positions.

```bash
sged-translate-coords \
    --sged HOG000003295_scores.sged \
    --index HOG000003295_SeqIndex.txt \
    --output HOG000003295_scores_ref.sged \
    --name AlnPos \
    --reverse
```

## Merge with evolutionary rate data

We use the translated DISEMBL file to merge it with site evolutionary rates, using the alignment coordinates as common index:

```bash
sged-merge \
    --sged1 HOG000003295_scores_ref.sged \
    --group1 AlnPos \
    --sged2 ../data/HOG000003295_charge_siteInfos.csv \
    --group2 Group \
    --output HOG000003295_scores+rates.sged
```

We can then compare the evolutionary rate (measured as the posterior substitution rate) and the intrinsic disorder:

```r
dat <- read.table("HOG000003295_scores+rates.sged", header = TRUE)
cor.test(~PR+HOTLOOP, dat, method = "kendall")
```

```
	Kendall's rank correlation tau

data:  PR and HOTLOOP
z = 5.6308, p-value = 1.794e-08
alternative hypothesis: true tau is not equal to 0
sample estimates:
      tau 
0.2604834 
```

Intrinsically disordered regions evolve faster!

