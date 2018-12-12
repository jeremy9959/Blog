---
type: posts
layout: single
excerpt: Binning and mapping reads in Ginkgo
tags: bioinformatics
---

## Making bins 

In general terms, ginkgo's approach is to construct bins with the property that an equal number of sequences map uniquely to points in the bin.  These bins are constructed
by first, fixing three parameters: 

- N, the expected number of reads per bin
- k, the length of the sequences produced by the sequencing machine (150, for example)
- align, the name of the alignment program (bwa or bowtie)

Next, within each chromosome with sequence C, consider in succession the sequences C[i:i+k] and map them to the reference genome using the chosen alignment program.  Keep only those reads that:

- have a quality score >= 30
- map uniquely to the reference genome (indicated by a 'NM:i:0' in the sam/bam file output by the chosen aligner)

The first bin then runs from 0 to the end of the Nth uniquely mapped sequence; the next from the beginning of the N+1-st uniquely mapped sequence to the end of sequence 2N, and so on.

In fact, ginkgo is a little bit cleverer than this:

- it locates the centromere of each chromosome by finding the longest string of 'N' in the sequence
- it divides the portions of the chromosome on either side of the centromere into bins
- since the length of these portions isn't a multiple of N, it distributes the remainder over the bins

Finally, ginkgo works out the GC content of each bin and saves that as well. This is done simply by counting the GgCc characters
as a fraction of the length of the sequence (excluding N characters).  This information is stored in a separate file and used later
for bias correction.

See:
- [genome/scripts](https://github.com/robertaboukhalil/ginkgo/tree/master/genomes/scripts) where the code for this is located
- [buildGenome](https://github.com/robertaboukhalil/ginkgo/blob/master/genomes/scripts/buildGenome), the particular script that 
controls the process.


## Mapping reads to bins in Ginkgo (compare ```bedtools coverage```)

Ginkgo maps the reads to bins by assigning the read to the bin where the read begins.  This is in contrast to the ```bedtools coverage```
command, which maps a read to a bin if it overlaps that bin.  In practice, the differences are slight. In one experiment,
90% of the bin counts computed by ```bedtools``` agreed exactly with the counts from ginkgo, and the remaining 10% of the counts differed
by less than 1%.

More specifically, we can illustrate this with some specific calculations.
We work with the following files to begin:
 - ```SCW-11_Bone-marrow.bed``` is the initial data file computed by running samToBed on a bam file from the alignment output.
 - ```variable_500000_76_bwa``` is the file computed by ginkgo identifying bins with 500kb unique reads of length 76 bp based on bwa alignment.  This file can be [downloaded](http://qb.cshl.edu/ginkgo/uploads/hg19.original.tar.gz) along with the many others from the ginkgo home site.
 
The small C++ program ```binUnsorted.cpp```  in the ginkgo distribution maps reads to the bins. It takes as arguments:
    1. the file of bins, 
    2. the length of that file, 
    3. the bed input file, 
    4. the name of the sample to be put as the first line of the output file, 
    5. and the name of the output file.

```
$ $GINKGO_PATH/scripts/binUnsorted variable_500000_76_bwa `wc -l < variable_500000_76_bwa` \ 
    SCW-11_Bone-marrow.bed SCW-11_Bone-marrow SCW-11_Bone-marrow_mapped
```


```python
!~/GitHub/ginkgo/scripts/binUnsorted variable_500000_76_bwa `wc -l < variable_500000_76_bwa` SCW-11_Bone-marrow.bed SCW-11_Bone-marrow SCW-11_Bone-marrow_mapped
```


```python
import pandas as pd
import numpy as np
import seaborn as sns


```


```python
binCounts = pd.read_table('SCW-11_Bone-marrow_mapped',names=['counts'],skiprows=1)
binCounts.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>counts</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3811</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1342</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2080</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1553</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2351</td>
    </tr>
  </tbody>
</table>
</div>



Excluding read counts greater than 5000 (of which there are 28), we obtain a nice poisson/negative binomial appearing distribution of counts.


```python
sns.distplot(binCounts[binCounts['counts']<5000],kde=False)
```


<figure>
<img src="/assets/images/output_6_1.png">
</figure>


Let's retry this using ```bedtools coverage```.  First, the bin file ```variable_500000_76_bwa``` needs to be put in bed format.
Currently it is only a list of the consecutive ends of the bins.  This can be done pretty easily with ```awk```.


```python
!cat variable_500000_76_bwa | awk 'BEGIN {x=0} ; /chr[0-9XY]+/ {if ($2>=x) { print($1"\t"x"\t"$2); x=$2 } else {x=0 ; print($1"\t"x"\t"$2); x=$2}}' > variable_500000_76_bwa.bed
```

Now we use bedtools coverage to compute counts.  We use the -F .75 flag so that a read must overlap a bin by at least 75% of its length to be counted.  Given the large bin sizes and the short reads this will make sure each read gets assigned to only one bin.  



```python
!bedtools coverage -counts -F .75 -a variable_500000_76_bwa.bed -b SCW-11_Bone-marrow.bed > SCW-11_Bone-marrow.counts
Counts = pd.read_table('SCW-11_Bone-marrow.counts',names=['chr','start','end','bedtools_count'])
```


```python
Counts['ginkgo_counts']=binCounts['counts']

```


```python
Counts.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chr</th>
      <th>start</th>
      <th>end</th>
      <th>bedtools_count</th>
      <th>ginkgo_counts</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>0</td>
      <td>978962</td>
      <td>3811</td>
      <td>3811</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>978962</td>
      <td>1483420</td>
      <td>1339</td>
      <td>1342</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>1483420</td>
      <td>1989002</td>
      <td>2080</td>
      <td>2080</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>1989002</td>
      <td>2493106</td>
      <td>1553</td>
      <td>1553</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1</td>
      <td>2493106</td>
      <td>3045688</td>
      <td>2351</td>
      <td>2351</td>
    </tr>
  </tbody>
</table>
</div>



The difference between the two counts is very slight -- less than 1%.  The two counts are identical 89% of the time, and never differ by even 1%. 


```python
Counts['diff']=Counts['bedtools_count']-Counts['ginkgo_counts']
print('Percent of bins where the counts differ:',100*(Counts['diff']!=0).sum()/Counts.shape[0],'percent')
print('Largest percent difference (in absolute value):',100* np.abs((Counts['diff']/Counts['bedtools_count'])).max(),"percent")
```

    Percent of bins where the counts differ: 11.509501613481534 percent
    Largest percent difference (in absolute value): 0.7905138339920948 percent






