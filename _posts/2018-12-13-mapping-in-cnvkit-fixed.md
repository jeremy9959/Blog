---
layout: single
title: 'Mapping to bins for cnvkit'
excerpt: 'cnvkit uses samtools bedcov for mapping'
tags: bioinformatics
---

```cnvkit``` uses, by default, the ```samtools bedcov``` command to compute the coverage.  This tool takes a bed file of regions and one (or more) bam files.  It returns, for each region (bin), the sum of the per-base read counts in that region.  The read depth computed by ```cnvkit``` is this sum divided by the length of the region.

By contrast, ```ginkgo``` counts the number of reads that start in each bin, and as we've seen this gives a number that is close to that returned by ```bedtools coverage```.   We would expect the ```ginkgo``` counts G, times the read length, divided by the interval length, to give the  ```cnvkit``` depth.


```python
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import linregress
```

The commands for comparison are:
```bash
$ cnvkit.py coverage $BAM_DIR/SCW-11_Bone-marrow.bam access-10kb.target.bed
$ bedtools coverage -F 1.0 -a access-10kb.target.bed -b SCW-11_Bone-marrow.bed > SCW-11_Bone-marrow.bedtools_coverage_F1.counts
```


```python
bedtCounts = pd.read_table('SCW-11_Bone-marrow.bedtools_coverage_F1.counts',header=None,names=['chr','start','end','-','count','x','y','z'])
bedtCounts = bedtCounts[['chr','start','end','count']]
cnvkCounts = pd.read_table('SCW-11_Bone-marrow.access-10kb.target.cnn')
```


```python
bedtCounts['depth'] = bedtCounts['count']/(bedtCounts['end']-bedtCounts['start'])
```


```python
linregress(bedtCounts['depth'],cnvkCounts['depth'])
```




LinregressResult(slope=70.32286807155764, intercept=0.00589153245627233, rvalue=0.9965766968261753, pvalue=0.0, stderr=0.07675388085434531)



This linear regression shows that the 'effective read length' is about 70 bp.

Just to show that in fact the cnvkit numbers agree with ```samtools bedcov``` let's compare.  The command is
```bash
$ samtools bedcov access-10kb.target.bed ~/hard_disk/SCG003/bam/SCW-11_Bone-marrow.bam > SCW-11_Bone-marrow.samtools_bedcov.counts
```

A quick linear regression shows that this agrees with the cnvkit coverage computation.


```python
samtCounts = pd.read_table('SCW-11_Bone-marrow.samtools_bedcov.counts',header=None,names=['chr','start','end','-','count'])
samtCounts['depth']=samtCounts['count']/(samtCounts['end']-samtCounts['start'])
linregress(samtCounts['depth'],cnvkCounts['depth'])
```




LinregressResult(slope=1.0000003802899826, intercept=-6.254817400130896e-08, rvalue=0.999999999999773, pvalue=0.0, stderr=8.865755503055448e-09)


