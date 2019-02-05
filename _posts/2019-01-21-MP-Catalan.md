---
type: posts
layout: single
excerpt: Some remarks on the moments of the Marchenko-Pastur distribution
title: Remarks on the Marchenko-Pastur Distribution
tags: mathematics
---

The proof that the Marchenko-Pastur distribution gives the asymptotic distribution of eigenvalues of
the Wishart Matrix of an $n\times p$ random matrix relies on combinatorial arguments similar to those
for the square case, where the moments are the Catalan numbers.

The argument in the MP case also relates the $k^{th}$ moment to the ordered trees on $k$ vertices.
However, it counts those trees differently.

The counting for the MP case first requires us to take the set of trees on $k$ vertices and mark the
nodes as "even" or "odd".  The root of the tree is always even, and the neighbors of an even vertex are always odd.
Having made such a labelling, we can partition the set of ordered trees into disjoint subsets depending on the number 
of even vertices.  

Among the 14 trees on 5 vertices, for example, the partition has the following distribution:

r|                                       1 | 2 | 3 | 4 |
--- | --- | --- | --- | ---|
number of trees with r even vertices | 1 |6 | 6 | 1 |


The tree with only $1$ even vertex is the star with 4 vertices meeting the root; the one with $4$ is the star
where the root is one of the leaves.

If $N(k,r)$ is the number of trees on $k$ vertices with $r$ even vertices, then the $k^{th}$ moment
of the Marchenko-Pastur distribution is

$$
\sum_{r=1}^{k-1} (p/n)^{r-1}N(k+1,r)
$$

where $p/n$ is the limiting number of rows/columns in the matrix.

As a check, consider the case where $\gamma=.5$.  Then the $4^{th}$ moment should be
$1+6(.5)+6(.5)^2+(.5)^3=5.625$.  A computation in sage confirms this.



```python
from sage.symbolic.integration.integral import definite_integral
R=RealField(200)
p=R(pi)
gamma1=(1+sqrt(.5))**2
gamma2=(1-sqrt(.5))**2
gamma=.5
var('x','f')
f=sqrt((x-gamma2)*(gamma1-x))/x/p/2/gamma
print('Check MP Formula: integral should be one (total mass): {}'.format(definite_integral(f,x,gamma2,gamma1)))
f4=x**4*f
print('4th moment is {}'.format(definite_integral(f4,x,gamma2,gamma1)))
```

    Check MP Formula: integral should be one (total mass): 0.9999999871272136
    4th moment is 5.624999941242612


A clever combinatorial argument gives a formula for these components--sort of "partial Catalan numbers".
The number of trees on $k$ vertices with $r$ even vertices $(1\le r\le k-1)$ is:

$$
N(k,r) = \frac{1}{r}\binom{k-1}{r-1}\binom{k-2}{r-1}
$$

A computation with $k=5$ and $r=1,2,3$ gives the answers above. For example
$$
N(5,3) = \frac{1}{3}\binom{4}{2}\binom{3}{2} = \frac{(1)(3)(6)}{3}=6.
$$


