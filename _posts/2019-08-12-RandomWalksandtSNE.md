---
layout:single
title: Random Walks on Graphs and tSNE
excerpt: a quick look at how tSNE uses random walks on graphs to compute affinities
---

## Random Walks and tSNE

One proposed technique to accelerate the tsne computation is to select a subset of the data points (landmark points) and work only with those.  It's important to capture the local neighborhood structure of the landmark points embedded in the full set of points, and the original tsne paper proposes a random walk method for doing this. First, one computes a neighborhood graph of the full data set. It's not clear to me how this is done, but leave it for now.  Then he affinities (the $p_{ij}$) are set to the probabilites that a random walk along the neighbor graph that leaves one landmark point first reaches another one.

### Harmonic functions and landmark points
Computing these probabilities is an exercise in graph theory.  The basic theorem is that if you want to find the probability that a random walk leaving a vertex $v_0$ first reaches a target vertex $x_0$ before reaching, say, $x_1,\ldots, x_n$,  you have to construct a harmonic function on the graph that takes the value one at $x_0$ and is zero at the other landmark vertices.


```python
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
G = nx.generators.random_graphs.gnp_random_graph(10,.3, seed=101)
nx.draw_networkx(G)
```


![png](/assets/images/RWtsne.png)


More generally, let's consider the "boundary value problem" on the graph, where we specify the values of our harmonic function at the "landmark points" and then try to find the function of minimum length with those specified values.  Decompose the linear span of the nodes into a "landmark subspace" $L$ and a regular subspace $N$.  
Bear in mind that the spaces $L$ and $N$ are not orthogonal with respect to $D$.

We have a fixed vector $\ell\in L$ and we want to find $n\in N$ so that
$\ell+n$ is of minimal length, in other words, we want to minimize
$$
H(n) = (\ell+n)^t D (\ell+n) = \ell^t D \ell + 2\ell^t D n + n^t D n
$$

Let's assume we are using the standard bases, so that $L$ corresponds to the first $k$ coordinates and $N$ to the remaining $N-k$ coordinates.  A harmonic function that is one at node $i$ and zero at the other landmark nodes corresponds to an $\ell$ that is $1$ in position $i$ and zero elsewhere.  The term $\ell^tDn$ is then the dot product of the last $N-k$ positions of the $i^th$ row of $D$ with the unknown vector $n$; call this $b_i$.  The term $n^t D n$ is the dot product computed using the lower $N-k \times N-k$ block of $D$; call this $D_*$.  The solution to this quadratic optimization problem is 

$$ f = -D_*^{-1}b_i.$$


To get the probabilities you must normalize 

For the graph above, here's an example, assuming we have nodes $0$ and $1$ as landmarks. First, find the Laplacian.



```python
D=nx.linalg.laplacian_matrix(G)
D=sp.sparse.csr_matrix.todense(D)
```

Next, find $-D_*^{-1}b_i$ for $i=0$:


```python
Q=-np.dot(np.linalg.inv(D[2:,2:]),D[2:,0])
Q
```

Notice that the points that are more likely to reach node zero before node 1 are points 2, 7, 8, and 9, which are the nodes close to 0.  

For node 1:


```python
Q=-np.dot(np.linalg.inv(D[2:,2:]),D[2:,1])
Q
```




    matrix([[0.1996419 ],
            [0.80662489],
            [0.61324978],
            [0.70993733],
            [0.68845121],
            [0.3992838 ],
            [0.37690242],
            [0.41987466]])



Here the more likely nodes are 3, 4, 5, 6 which are the ones close to node 1.

## tSNE complications

The only flaw in the discussion above is that, for tSNE purposes, we want to consider paths that not only END on a landmark point but also begin on one.  The paper proposes duplicating the landmark nodes and considering paths that start on one of the duplicates.  

It's not clear to me what this actually means, so we'll have to come back to it.
