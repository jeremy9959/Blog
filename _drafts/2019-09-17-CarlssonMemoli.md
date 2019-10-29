---
title: Carlsson-Memoli on Hierarchical Clustering
type: posts
layout: single
excerpt: Carlsson and Memoli show that single-linkage hierarchical clustering has good properties (and other types do not)
tags: clustering
---

The paper [Characterization, Stability, and Convergence of Hierarchical Clustering Methods](http://www.jmlr.org/papers/volume11/carlsson10a/carlsson10a.pdf)
by Carlsson and Memoli applies topological ideas to study hierarchical clustering.  They obtain a result complementary to that of Kleinberg,
in the sense that they show single-linkage hierarchical clustering does have good categorical properties as a functor from finite metric spaces to ultrametric trees.

By a partition of a finite set $X$ we mean a disjoint covering of $X$ by subsets.

**Defintion:** Let $X$ be a finite set and $\theta:[0,\infty]\to \mathcal{P}(X)$ where $\mathcal{P}(X)$ is the set of partitions of $X$.  The pair $(X,\theta)$
is a dendrogram if:

- $\theta(0)$ is the maximal partition of $X$ into one-element subsets.
- For sufficiently large $t$, $\theta(t)$ is the trivial partition of $X$ into a single subset.
- If $r<s$, then $\theta(r)$ is a refinement of $\theta(s)$.
- For all $r$ there exists $\epsilon>0$ so that $theta(r)=\theta(s)$ for all $s\in [r,r+\epsilon]$

A dendrogram on $X$ is equivalent to an ultrametric on $X$; that is, a metric that satisfies the stronger triangle inequality $\rho(x,y)\le\max(\rho(x,z),\rho(z,y))$.
An ultrametric is equivalent to a dendrogram; given an ultrametric, the associated dendrogram comes from the partitions whose sets
are the equivalence classes $\rho(x,y)\le r$.  Given a dendrogram, the associated ultrametric is setting $\rho(x,y)$ to be the smallest value of $r$ for which
the two points $x$ and $y$ are in the same element of the partition.

A hierarchical clustering algorithm associates a dendrogram to a finite metric space.  In fact, Carlsson-Memoli's
definition of a hierarchical clustering method is a map The simplest (and in some sense, as C-M show, the most canonical)
such algorithm is "single linkage clustering" which has several definitions.  


**Definition 1:** For each $r$, define an equivalence relation on $X$ so that $x\sim_{r}y$ when there is a sequence of points $x_0=x,x_1,\ldots,x_k=y$
so that all of the  distances $\rho(x_i,x_{i+1})$ are less than or equal to $r$.  Equivalently, construct an ultrametric by settings
$$
\mu(x,y)=\inf\max\{\rho(x_i,x_{i+1})\}
$$
where the maximum is taken over all the steps in a sequence of points running from $x$ to $y$, and the infimum is taken over all such sequences;
then take the associated dendrogram.  This is called the **maximum subdominant ultrametric.** It is the largest among all ultrametrics that are
less than or equal to the original metric $\rho$. 

One feature of this construction which may be non-standard is that it allows multiple points to coalesce at the same time.  For example,
consider a finite graph where all the edge lengths are one and the distance is the length of the shortest path between vertices.  The associated
ultrametric puts all points at distance one from each other, so the associated dendrogram collapses everything at time one. *This is different than
the single linkage algorithm in the ```sklearn``` library, for example, which constructs a binary tree by making arbitrary choices of what to merge
when there is a group of equidistant points.*


**Theorem: (C-M)**  Let $\mathcal{T}$ be a hierarchical clustering method -- that is, a map which associates to a finite metric space $X$
an ultrametric space with the same underlying points.  Then the following properties characterize single linkage clustering:

- On the two point set, where $x$ and $y$ are at distance $\delta$, one obtains the two point set with the ultrametric $u$ such that $u(x,y)=\delta$.
- Whenever $\phi X:\to Y$ is a distance non-increasing map, the induced map on ultrametrics is also distance non-decreasing.
- The mimimum non-trivial distance between points of $X$ measured by the ultrametric is at least the minimum non-trivial distance measured by the metric.

Implicit in this paper by C-M is a functorial notion of hierarchical clustering that is made more explicit in their
paper [Classifying Clustering Schemes](https://link.springer.com/article/10.1007/s10208-012-9141-9).  
Let  $\mathcal{M}$ be the category of finite metric spaces with distance non-increasing maps, and let $\mathcal{P}$ be the category
of "persistent sets" -- where a persistent set is defined exactly like a dendrogram but without the condition that there is a large enough $t$
such that $\theta(t)$ is all of $X$.  This means that one has an ultrametric but some of the distances are infinite. 

Then the rule that sends a finite metric space $X$ to the same set with the maximum subdominant ultrametric is a functor from the
category of finite metric spaces with distance non-decreasing maps to the category of "persistent sets", where a morphism $f:X\to Y$
has the property that $f^{*}(\theta_{Y}(r))$ is a refinement of $\theta_X(r)$ for all $r$. Essentially this means that the map is distance
non-increasing between the metric and the ultrametric.

What is particularly interesting about this is that other methods are NOT functorial.  For example, one can take what's called complete linkage,
or average linkage
ein which the distance between clusters is given by the *largest* distance between a pair of points, one in each cluster, or the average distance
between points in the two clusters. Carlsson-Memoli
show that this rule is NOT functorial, which in practice means that slight perturbations of the metric yield very different dendrograms.  Thus
single linkage has good technical properties EVEN THOUGH in practice it can take a long string of points in a line and put them in a single
cluster.

The last part of the C-M paper considers maps between different metric spaces, and how that relates to the clustering.  They use the *Gromov-Hausdorff*
distance to compare the original metric spaces and the associated ultrametric spaces.  They prove that the functor giving single linkage
clustering reduces the GH distance.  The most interesting part of this is a beautiful result.

**Theorem:** (C-M, Theorem 28) Let $Z$ be a compact metric space.  Let $X$ and $X'$ be any two finite subsets of $Z$ and let $Y$ and $Y'$ be
the associated ultrametric spaces obtained by the single-linkage clustering functor.  Then:

- The GH distance between $Y$ and $Y'$ is bounded above by $d_{H}^{Z}(X,Z)+d_{H}^{Z}(X',Z)$ where $d_{H}^{Z}$ is the Hausdorff distance in $Z$.
- Suppose that $Z$ decomposes into compact, disjoint, path-connected components $Z_{1},\ldots, Z_{m}$.  Let $A$ be the finite metric space whose points
are the components $Z_{i}$ of $Z$ and whose distances are by the Hausdorff distance between compact sets. Then if $d_{H}^{Z}(X,Z)<\mathrm{sep}(A,d_{A})/2$,
we know that the GH-distance between $Y$ and $A$ with its ultrametric is at most $d_{H}^{Z}(X,Z)$. 
- If $X_{n}$ is a sequence of finite subsets of $Z$, with the induced metric, and so that $d_{H}^{Z}(X_n,Z)\to 0$ as $n\to\infty$, Then the GH distance
between $X_n$ with the single linkage ultrametric and $A$ with the ultrametric goes to zero.

In other words, the clustering structure of finite subsets of a compact space captures, in the limit, the clustering structure of the components.

Question: does the fact that this all takes place in a compact space avoid the chaining phenomenon?




