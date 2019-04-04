---
type: posts
layout: single
title: Samtools and Pysam
excerpt: Reconciling output from samtools mpileup and pysam mpileup
tags: bioinformatics python
---

I struggled quite a bit to reconcile the output from pysam and from samtools mpileup.
This little piece of code reads a sorted bam file using the pileup API for pysam
and also runs samtools mpileup and does a comparison
The key subtleties seem to be:
- the option for BAQ filtering is on by default in the samtools mpileup 
- one has to parse the samtools sequences carefully to handle insertions, deletions
- its import to use the 'indel' flag in pysam 
- pysam uses zero based coordinates while samtools uses one-based
- samtools can set the depth but can pysam?
- In the sequence string output by samtools, a ^ means the start of a read sequence, but the
- next character is a base quality -- NOTE that it could happen to be an A, for example, but it's not - a base, it's a quality, and you shouldn't count it.


First load the libraries.

```python
import pysam
import pandas as pd
from collections import Counter
import subprocess
import re
```

There is nothing special about this file, it's just useful for testing. 

```python
sampath = '/home/jet08013/usb_disk/cell_paper_bam/Processed/SRR8193067.sorted.bam'
samfile = pysam.AlignmentFile(sampath,'rb')
```

For these options:

```stepper = nofilter``` means read everything

```truncate = False``` (the default) means that you count all reads that overlap the target region
we turn this off, though here there's no region specified so it doesn't matter.

```python
MTpileup = samfile.pileup('MT',stepper='nofilter',max_depth=500000,truncate=False,min_base_quality=0)
df_pysam = pd.DataFrame(np.zeros(shape=(16569,5)),columns=['C','A','G','T','isindel'],index=list(range(16569)))
df_samtools = pd.DataFrame(np.zeros(shape=(16569,5)),columns=['C','A','G','T','isindel'],index=list(range(16569)))
```

First, we do the pysam approach

```python
for position in MTpileup:
    c=Counter({'C':0,'A':0,'G':0,'T':0,'isindel':0}) 
    for pileupread in position.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            c[pileupread.alignment.query_sequence[pileupread.query_position].upper()]+=1
        if pileupread.indel:
            c['isindel']+=1
    df_pysam.loc[position.reference_pos+1]=pd.Series(c)  # note that we add 1 to compare with samtools
samfile.close()
```

Next, we do the samtools subprocess approach.  We send ```stdout=subprocess.PIPE``` so we can read the output from the file,
but we must have ```universal_newlines=True``` for this to work correctly.

```python
subproc=subprocess.Popen(['samtools', 'mpileup','-d','500000','/home/jet08013/usb_disk/cell_paper_bam/Processed/SRR8193067.sorted.bam'],stdout=subprocess.PIPE,universal_newlines=True)
```

No reference is specified, so the sequence contains the actual base pairs, annotated with '^' and '$'.

```python

for x in subproc.stdout:
    c=Counter({'C':0,'A':0,'G':0,'T':0,'isindel':0})
    (chrom,loc,ref,depth,sequence,quals)=x.split() 
    i=0
    while i<len(sequence):
        if sequence[i] not in ['+','-', '^']:
            c[sequence[i].upper()]+=1
            i+=1
            continue
        else:
            if sequence[i]=='^':
                i+=2
            else:
                i+=1
                indel_len=re.match('[0-9]+',sequence[i:])[0]
                c['isindel']+=1
                i+=int(indel_len)+len(indel_len)
    c[ref]=c[',']+c['.']
    df_samtools.loc[int(loc)]=pd.Series(c)
```

Check to see if we are on target:


```python
for base in ['C','A','G','T','isindel']:
    print(base, (df_samtools[base]!=df_pysam[base]).sum())
```

We should get all zeros, and it seems that we did!

