# sgedtools
Tools to handle files in the SGED (Site Group Extensible Data) format

The SGED format allows to annotate (group of) sites in a sequence alignment.
It is essentially a tabular format (TSV or CSV), with one group per line.
The first column describes the site group, as 1-based positions separated by a semi-colon (';'), within square brackets.
Other columns contains several annotation values.

Example, as output by the program CoMap:
```
Group	Size	IsConstant	Dmax	Stat	Nmin	p.value	nobs	code	Method	FDR
[147;157]	2	no	0.183705	0.816295	1.56728	0.0437622955991075	348701	*	charge	no
[353;397]	2	no	0.443686	0.556314	3.23407	0.0461919112116152	211378	*	charge	no
[367;413]	2	no	0.122811	0.877189	1.00574	0.0479694948963177	273527	*	charge	no
[143;146]	2	no	0.526655	0.473345	6.43379	0.0533397809221454	47562	.	charge	no
```

This package contains a series of scripts to manipulate such file.
These scripts are written in Python, and make (sometimes) use of the BioPython libraries (for sequence alignment and 3D structure manipulation).
