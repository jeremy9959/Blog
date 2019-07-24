---
type: post
layout: single
title: Counting Haplotypes using Expectation Maximzation
excerpt: notes on the 1995 paper by Excoffier and Slatkin
tags: math bioinformatics
---

## Counting Haplotypes

Some notes on the classic paper

[Excoffier, M. and Slatkin, L. Maximum-Likelihood Estimation of Molecular Haplotype Frequencies in a Diploid Population. Mol. Biol. Evol. 12(5):921-927 (1995)](http://www.uvm.edu/~rsingle/stat380/F04/possible/Excoffier+slatkin-MolBiolEvol-1995.pdf)

Suppose we have two loci, each of which has potential genotype at that locus of  $0,1,2$ representing the number of "reference" alleles.  An individual can have one of nine combined genotypes depending on the situation at each locus.

One can ask a more refined question: what are the underlying haplotypes?  If the individual is, say, of type $22$ then it must be the case that each chromosome is of type $11$.  But if the individual is of type $11$
then there are two possibilities: $(10,01)$ or $(00,11)$.

This question is of interest because one might want to know if the sites are independent of one another or if there is non-trivial linkage between them.

More generally, if there are $n$ loci in a diploid individual, there are $3^n$ possible combinations of genotypes.  Although it is somewhat of an abuse of language, we follow the language of the paper above and call such a combination a phenotype.  Phenotypes are indexed by vectors of $n$ integers between $0$ and $2$.

We have $2^n$ haplotypes obtained by choosing either the reference or alternate allele at each of the n sites,
so haplotypes are indexed by vectors of $n$ integers between $0$ and $1$.  The phenotype associated with a pair
$\alpha$ and $\beta$ of haplotypes corresponds to the vector sum $\alpha+\beta$.  If a phenotype has a zero or a two at a particular location, then that determines the possible haplotypes at that location.  If the phenotype has a $1$, then at that location there are two possibilities.  As a result, the number of ways that a phenotype
can be decomposed into a sum of haplotypes (taking symmetry into account) is $2^{s-1}$ where $s$ is the number of $1$'s in the phenotype (or $1$ if that number is zero).

### EM algorithm

We're given counts of each of the $3^n$ phenotypes and we want to give the maximum likelihood estimate of the haplotype frequencies under the assumption that there is no linkage between them.

Let's  consider the case where $n=2$ so that there are $9$ phenotypes and $4$ haplotypes.   The data consists of counts $N_{ij}$ over the nine phenotypes (so $i$ and $j$ are between $0$ and $2$ inclusive). Let $N$ be the sum of the $N_{ij}$.

We choose initial (prior) estimates for the haplotype frequencies $p_{ij}^{(0)}$ (where here $i$ and $j$ are between $0$ and $1$.)

In light of the phenotype data, we find the expected number of times each haplotype occurs assuming the
haplotype frequency estimates.  In most cases, the phenotype determines the haplotypes, so that
for example phenotype $22$ contributes two copies of haplotype $11$ and phenotype $12$ contributes
$1$ copy of haplotype $01$ and one copy of haplotype $11$.  The only exception is phenotype $11$, which
is made up of either haplotype $00$ and haplotype $11$ or haplotype $10$ and haplotype $01$.  So among
the observed number of individuals with phenotype $11$, how many come from each haplotype?  Given
the estimated haplotype frequencies, we can compute the expected fractions as
$$
E(10,01)=\frac{p_{10}p_{01}}{p_{01}p_{10}+p_{00}p_{11}}
$$
and
$$
E(00,11)=\frac{p_{00}p_{11}}{p_{01}p_{10}+p_{00}p_{11}}.
$$

Putting this information together we can compute the estimated number of times each haplotype appears given the data.  Dividing by the total number of chromosomes, we get the expected frequency of each haplotype assuming
the data and the prior estimates.

Now we view these frequencies as the maximum likelihood estimate of the true probabilities given the "expected" data.  Now we put those frequencies back through the process to refine them further.






### Numerical Example

Let's suppose that we observe the following number of phenotypes where entry $i,j$ is the count of phenotype
$ij$.

|  N| 0 | 1 | 2 |
|---|---|---|---|
| 0 | 3 | 17| 17|
| 1 | 51|158|98 |
| 2 | 97|340|219|

Our unknown parameters are the haplotype frequencies $p_{00},p_{01},p_{10},p_{11}$.

In every case except for the case of $N_{11}$, we know that the phenotype determines the genotype.
But the $N_{11}$ case is a combination of genotypes $(00,11)$ and $(10,01)$.

Let's suppose our initial guess for the $p_{ij}$ is that all of them are equal to $.25$. Under the assumption that an individual genotype is obtained by choosing the haplotypes independently at random, then among the $158$
individuals with genotype $11$ we would expect $79$ of each type.

The maximum likelihood estimate for the $p_{ij}$ then comes by counting the observed frequency of each haplotype.
For example, consider the $00$ haplotype.  The genotype $00$ contributes two such haplotypes.  The genotypes $01$ and $10$ each contribute one such haplotype.  Among the $158$ examples of $11$ genotypes, our expectation is that 79 are of type $(00,11)$ and these contribute 1; the other possibility contributes none. The other genotypes don't include any haplotypes of type $00$.  This means that
under our proposal, we see $2\times 3 + 17 + 51 + 79=153$ haplotypes of type $00$ out of a total number of $2000$ observed. (The sum of the entries in the matrix above is $1000$).  Thus the new estimate for $p_{00}=.153/2$.

We repeat this process for each of the other three $p_{ij}$ to complete the $M$-step; and then we recalculate the
distribution of the $158$ entries of type $N_{11}$ among the possible genotypes.


To construct a complete numerical example, let's choose haplotype probabilities and compute a sample matrix N.



```python
import numpy as np
from numpy.random import multinomial
np.random.seed(8) # for reproducibility
(p00,p01,p10,p11)=.05,.13,.35,.47
P=[p00**2,2*p00*p01,p01**2,2*p10*p00,2*p00*p11+2*p10*p01,2*p01*p11,p10**2,2*p10*p11,p11**2]
N=multinomial(1000,P).reshape((3,3))
print('Observed phenotype counts\n',N)
```

Observed phenotype counts
[[  4  20  21]
[ 31 128 125]
[127 324 220]]



```python
def EM(p00,p01,p10,p11):
n00=2*N[0,0]+N[0,1]+N[1,0]+p00*p11/(p00*p11+p10*p10)*N[1,1]
n10=N[1,0]+p10*p01/(p00*p11+p10*p01)*N[1,1]+N[2,1]+2*N[2,0]
n01=N[0,1]+p10*p01/(p00*p11+p10*p01)*N[1,1]+N[1,2]+2*N[0,2]
n11=2*N[2,2]+p00*p11/(p00*p11+p10*p10)*N[1,1]+N[1,2]+N[2,1]
return n00/2000,n01/2000,n10/2000,n11/2000
```

Iterating from our initial guess we obtain the sequence of improvements.


```python
guess=(.25,.25,.25,.25)
for i in range(10):
print('{:.3f} {:.3f} {:.3f} {:.3f}'.format(*guess))
guess=EM(*guess)
```

0.250 0.250 0.250 0.250
0.061 0.126 0.337 0.476
0.043 0.131 0.342 0.458
0.039 0.138 0.349 0.454
0.038 0.140 0.351 0.453
0.037 0.141 0.352 0.452
0.037 0.141 0.352 0.452
0.037 0.141 0.352 0.452
0.037 0.141 0.352 0.452
0.037 0.141 0.352 0.452


Notice that the last line is close to our original values of .05,.13,.35,.47.
