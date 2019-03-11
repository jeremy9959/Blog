---
type: posts
layout: single
tags: bioinformatics
title: Computing Autocorrelation
excerpt: Computing autocorrelation of long sequences
---

In looking at the work of Zhang, et. al. in [this post](/AutoCorrelation) and [this one](/AutoCorrelation2), I ran into the problem
of computing the autocorrelation of a very long sequence.

For example, the coverage data used in Figure 2a reproduced [here](/AutoCorrelation2) is 4GB of lines that look like this:
```
woodmont $ head SRR1779331.depth
CM000663.1	12193	1
CM000663.1	12194	1
CM000663.1	12195	1
CM000663.1	12196	1
CM000663.1	12197	1
CM000663.1	12198	1
CM000663.1	12199	1
CM000663.1	12200	1
CM000663.1	12201	1
CM000663.1	12202	1
...
```

One can read this into a pandas dataframe, and save memory by omitting the first column (since it just denotes the first chromosome and it's the same on every line) and treating the second (the location) and the third (the depth) columns as integers.  Nevertheless
the resulting dataframe is 2.8GB and takes a while to load.

After experimenting with a few ideas, the simplest memory-efficient approach seems to be to use the sparse matrix routines
in scipy to compute the autocorrelation.  Here is the code:

```python
import scipy.sparse as sp
import pandas as pd
import numpy as np
from tqdm import tqdm

df = pd.read_csv('/home/jet08013/hard_disk/ZhangPaper/depths/SRR1779331.depth',sep='\t',usecols=[1,2],header=None,names=['loc','depth'])
```

We need to know the length of chromosome 1 in order to compute means and so on.

```
CHR1_LENGTH = 249250621
```

The are a number of constructors for a csc matrix, but
the one I use takes an array of data and then two arrays ```indices``` and ```indptr``` that specify where to 
put the data into the matrix.  The ```indices``` gives you the row position of the corresponding data element.  
The  ```indptr``` array allows you to figure out which column to put the data in.  Formally 
speaking, ```indices[indptr[i]:indptr[i+1]]```
holds the row indices for column i, and ```data[indptr[i]:indptr[i+1]]``` holds
the corresponding data.

An example might help. Suppose:

```
data = [1,3,4,6]
indices = [0, 0, 1, 2]
indptr=[0,1,1,4]
```

- ```[indptr[0]:indptr[1]]=[0:1] = [0]```, so ```data[0]=1``` goes in column zero.  
- Since ```indices[0]=0```, ```data[0]``` goes in position ```(0,0)```.
-  ```[indptr[1]:indptr[2]]=[2:2]``` which is empty.  So column 1 is zero.
- ```[indptr[2]:indptr[3]]=[1:4] = [1,2,3]```.  So ```data[1:3]``` goes in column 2, and since ```indices[1:3]=0,1,2```, these data points
go in position ```(0,2), (1,2), and (2,2)```.

The final matrix is ```[[1,0,3],[0,0,4],[0,0,6]].```

Of course we are in the much simpler case of a single column.

```
data = df['depth'].values
indices = df['loc'].values
indptr = [0, len(indices)]

M = sp.csc_matrix((data,indices,indptr),shape=(CHR1_LENGTH,1))
mu = M.sum()/M.shape[0]
```

The great thing about this representation is that you can "shift" the matrix for purposes of computing the autocorrelation
just by subtracting from the index, rather than actually  moving any data. 

Values that are shifted into negative indices are effectively set to zero, since the sparse routines ignore such entries.

The loop below computes the auto correlation
using Zhang's normalization for shifts of length $2^n$ and saves the results in the matrix acor.

```
acor = np.zeros(25)
for i in tqdm(range(0,25)):
    M1 = sp.csc_matrix((data,indices-(2**i),indptr),shape=(CHR1_LENGTH, 1))
    N = M.multiply(M1)
    acor[i]=(N.sum()/CHR1_LENGTH-mu**2)/mu**2
```

Experiments show that the memory usage of this code is constant through the loop.

## Roads not taken

There is a much more efficient approach that would compute the autocorrelation *as you read the file*. This would
require buffering only within the range of the longest lags under consideration, which in this case is on the order
of 1M locations (rather than 250M).  But the approach here, exploiting the fast numpy code, saves a lot of development time.



