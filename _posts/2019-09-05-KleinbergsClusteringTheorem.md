---
title: Kleinberg's Impossibility Theorem for Clustering
type: posts
layout: single
excerpt: Kleinberg shows that it's impossible to find a clustering algorithm that satisfies three simple properties.
tags: clustering
---

In the paper [An Impossibility Theorem for Clustering](https://www.cs.cornell.edu/home/kleinber/nips15.pdf),
Jon Kleinberg introduces three simple properties that one might hope a clustering algorithm would satisfy,
and then proves that no algorithm can satisfy all three.

Suppose given a set $S$ with $n\ge 2$ points and a distance function $d: S\times S\to \mathbf{R}$ such that $d(i,j)=0$ only if $i=j$
and such that $d(i,j)=d(j,i)$ for all $(i,j)\in S\times S$.  Note that we don't assume that, for example, $d$ is a metric,
but the theorem is true even if we restrict to that class of distance functions.

Let $\mathcal{D}(S)$ denote the set of distance functions on $S$ and let $\Pi(S)$ be the set of partitions of $S$ into disjoint subsets. 

**Definition:** A clustering function is a function $f:\mathcal{D}(S)\to \Pi(S)$; given a distance function, $f$ returns
a partition of $S$ into disjoint clusters.

Kleinberg considers three properties one might expect of a clustering function.

- *Scale Invariance*: This asserts that $f(d)=f(\alpha d)$ for all distance functions $d\in\mathcal{D}(S)$ and all real $\alpha>0$.
- *Richness*: Given a partition $\Gamma\in \Pi(S)$, there is a $d\in \mathcal{D}(S)$ so that $f(d)=\Gamma$. In other words, $f$ is surjective.
- *Consistency*: Suppose $d$ and $d'$ are two distance functions and let $\Gamma$ be partition of $S$.  We say that $d'$ is a $\Gamma$-transformation of $d$ if $d'(i,j)\le d(i,j)$ for all pairs $(i,j)$ belonging to the same cluster in $\Gamma$, while $d'(i,j)\ge d(i,j)$
for all pairs belonging to different clusters.  Then $d$ and $d'$ are consistent if, whenever $d'$ is an $f(d)$-transformation of $d$,
then $f(d')=f(d)$.

There are clustering functions that satisfy any two of the three
conditions.  For concreteness assume that $S$ are the nodes of a
graph, connected by edges of weight $d(i,j)$.  The clustering
functions find subgraphs of this graph by choosing a subset of the
edges according to a rule.

1.  Fix $1<k<n$. Put the edges in order by non-decreasing weight and add edges to the subgraph until it has exactly $k$-connected
components. Use lexicographic order to break ties. These components are the clusters. (This is agglomerative clustering)
2. Fix a distance $r$ and add all edges of weight at most $r$.  The connected components are the clusters.
3. Fix $1>\alpha>0$ and let $R$ be the maximum value of $d$.  Add edges of weight at most $\alpha d$.

**Proposition:** Method 1 satisfies Scale-invariance and Consistency; Method 2 satisfies Scale-Invariance and Richness;
Method 3 satisfies Richness and Consistency.

In case *1*, scaling the lengths doesn't affect their order so the same components are constructed.  To see consistency,
assume $\Gamma$ is the partition arising from $d$.  If $d'$ a $\Gamma$-transformation of $d$, it means that $d'(i,j)\le d(i,j)$
for all edges $i,j$ added to the subgraph, while $d'(i,j)\ge d(i,j)$ for all edges not yet in the subgraph.  Since we used
lexicographic order to break ties, and that doesn't depend on $d$ or $d'$, the two distance functions yield the same ordered list of edges and thus the same subgraph and the same clusters.  Since you never get more than $k<n$ clusters, you don't have richness.

In case *2*, you don't have scale invariance because changing $r$ changes the clusters.  If you want to get a particular cluster,
you can have all edges within a cluster have weight smaller than one, and all edges that cross clusters have weight greater than 1.
For consistency, suppose that $d$ gives rise to a particular partition and that $d'(i,j)\le d(i,j)$ within clusters and
$d'(i,j)\ge d(i,j)$ between clusters.  Then you end up adding exactly the same edges to the subgraph for $d'$, and therefore you get
the same $\Gamma$. 

In case *3*, it's clearly scale invariant and for a particular set $S$ you can choose a $d$ that is $1$ within the desired
clusters and greater than $1/\alpha$ between them.  Thus you can get any set of clusters you want.  To see that 
consistency fails, we need at least three points and so we have at least three edges.  Suppose that $d$ is the constant function
taking the value $1$.  Then the associated clusters are the distinct points.  Now choose just one pair of points and construct $d'$ with
the value $1/\alpha$ there while $d'=1$ everywhere else.  Now between clusters we will have $d'\ge d$.  However, the maximum is now $1/\alpha$ so the clustering algorithm joins all points at distance less than or equal to $1$; so we end up with all the points, except possible one of them, joined together.


**Theorem: (Kleinberg)** For each $n\ge 2$ there is no clustering function $f$ satisfying Scale-Invariance, Richness, and Consistency.

The proof follows from this theorem.

**Theorem:** If a clustering function satisfies scale-invariance and consistency, then the range of $f$ is an antichain -- meaning
that there is no pair $d$, $d'$ of distance functions so that $f(d)$ is a refinement  of $f(d')$.  Put another way, if $f(d)$ is
a refinement of $f(d')$, then $f(d)=f(d')$.  

To prove this second theorem, suppose $f$ is a clustering function that satisfies consistency and scale invariance.  Let $d'$
be a distance function
let $\Gamma'=f(d')$, let $d$ be another distance function, and suppose that $\Gamma=f(d)$ is a refinement of $\Gamma'$.

Let $a'$ be the minimum distance among points in the same cluster of $\Gamma'$, and let $b'$ be the maximum distance
among points in different clusters of $\Gamma'$.  Choose $a$ and $b$ similarly for $\Gamma$. Consistency tells us
that if $d^{\dagger}$ is a distance function that is less than $a'$ within clusters of $\Gamma'$, and greater than $b'$ between them,
then $d^{\dagger}=d'$; and similarly for $a$, $b$, $\Gamma$, and $d$. 

Since $\Gamma$ is a refinement of $\Gamma'$, if $i,j$ are in different clusters of $\Gamma'$ they must be in different clusters of 
$\Gamma$.  

Choose $0<\epsilon<aa'/b$.  Define a new distance function $d^{\dagger}$ with the following properties:

- $d^{\dagger}(i,j)=\epsilon$ if $i,j$ are in the same cluster of $\Gamma$
- $d^{\dagger}(i,j) = a'$ if $i,j$ are in the same cluster of $\Gamma'$, but not of $\Gamma$.
- $d^{\dagger}(i,j) = b'$ if $i,j$ are in different clusters of $\Gamma'$. 

The second two properties say that $d^{\dagger}$ is consistent with $d'$, so $f(d^{\dagger})=\Gamma'$.  On the other hand, let $\alpha = b/a'$.
Then $\alpha d^{\dagger}$ satisfies

- $\alpha d^{\dagger}(i,j) = (b/a')\epsilon < a$ if $i,j$ are in the same cluster of $\Gamma$
- $\alpha d^{\dagger}(i,j) = (b/a')a' = b$ if $i,j$ are in the same cluster of $\Gamma'$ but not of $\Gamma$.
- $\alpha d^{\dagger}(i,j) = (b/a')b' \ge b$ since $b'\ge a'$ if $i,j$ are in different clusters of $\Gamma'$.

The upshot of this is that $\alpha d^{\dagger}$ is gives the same partition as $d$ by consistency, and the same partition as $d'$
by scale invariance -- in other words, $\Gamma=\Gamma'$. 

This proof does not try to preserve any additional properties of $d$, such as the triangle inequality; Kleinberg shows that
one can improve the inequalities and preserve such additional conditions without affecting the result.

To fully characterize possible partitions, Kleinberg proves that any antichain can arise via a clustering algorithm
satisfying all three conditions.

**Theorem:** Given any antichain of partitions (that is, any set of partitions, none of which are refinements of another one),
there is a clustering function whose range is that antichain.

Fix an antichain $\mathcal{A}$ and consider the objective function
$$
\Phi_{d}(\Gamma) = \sum_{(i,j)\sim\Gamma} d(i,j)
$$
where $(i,j)\sim\Gamma$ means $i$ and $j$ are in the same subset of $\Gamma$.  Let $\Gamma(d)$ be the partition **in** $\mathcal{A}$
that minimizes this objective function.  Kleinberg shows that this clustering function is scale-invariant and consistent,
and has range $\mathcal{A}$.  The proof shows how to define $d$ for a given partition $\Gamma$, and then shows that
the result is consistent.   Note (as does Kleinberg) that the minimization must consider only partitions in $\mathcal{A}$ or
one will obtain the trivial partition.







