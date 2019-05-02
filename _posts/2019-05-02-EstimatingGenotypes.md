---
type: posts
layout: single
title: Estimating genotypes and testing HWE from sequences
tags: bioinformatics
excerpt: Continuing discussion of the Li paper Inference using sequencing data...
---

## Estimating genotypes and testing for HWE

We continue to read the old paper 

Li, H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, Vol. 27, no. 21, 2011, pp. 2987-2993.

In the [earlier discussion](https://jeremy9959.github.io/genotype-likelihoods) we used the sequencing results to estimate the population frequency of the reference allele assuming that the population is in Hardy-Weinberg equilibrium.

Now we consider the problem of estimating the genotype frequencies independent of the HWE hypothesis.  Combining the two results then gives us a likelihood ratio test for HWE.

We recall the computation from the [earlier discussion](https://jeremy9959.github.io/genotype-likelihoods) which yielded the formula
\begin{eqnarray}
\mathcal{L}(g) & = & P(\hbox{$l$ reference in $N$ reads}|g) \cr
&=& \frac{1}{2^{N}}\prod_{j=1}^{l}\left((m-g)\epsilon_{j}+g(1-\epsilon_{j})\right)\prod_{j=l+1}^{N}\left((m-g)(1-\epsilon_{j})+g\epsilon_{j})\right) \cr
\end{eqnarray}


where $N$ is the depth of the sequencing, $\epsilon_{j}$ is the error probability of the $j^{th}$ read, and we assume ploidy = 2.

If $\overline{\xi}=(\xi_{0},\xi_{1},\xi_{2})$ are the frequencies of each genotype $g=0,1,2$ (which is what we are trying to estimate),
we have, for $S$ samples, that

$$
\mathcal{L}((\xi_{0},\xi_{1},\xi_{2})) = Pr(\hbox{$l$ reference in $N$ reads}|\overline{\xi}) = \prod_{i=1}^{S}\sum_{g=0}^{2}\mathcal{L}_i(g)\xi_{g}
$$

Note that $\sum_{i=0}^{2}\xi_{i}=1$ by definition.

The EM algorithm for this situation means that we choose an initial guess for the $\xi_g$.
Strictly speaking what we do is would put a delta prior on (for example) $\xi_0$ with $\xi_0=\xi_0^{*}$.  Then we compute
$$
E[P(\xi_0|d)]=\frac{1}{S}\sum_{i=1}^{S}\frac{P(d|\xi_0^{*})}{\sum_{g=0}^{2}P(d|\xi_g^{*})\xi_g^{*}}
$$
and then our next estimate for $\xi_i$ is to set $\xi_i=E[P(\xi_i|d)]$ and iterate this to obtain the maximum likelihood estimate for the genotypes.

To test for HWE, we compute a test statistic 
$$
D=-2\log \frac{L_\psi}{L_{\xi}}
$$
where $L_{\psi}$ is the log-likelihood obtained by setting 
$$(\xi_0,\xi_1,\xi_2)=((1-\psi)^2,2\psi(1-\psi),\psi^2)$$
using the maximum likelihood estimate for $\psi$ computed in the earlier discussion that assumes HWE.

The idea here is that the measure $D$ will be small if the genotype frequences are close to those predicted by the HWE hypothesis.  

More precisely, an asymptotic result says that $D$ follows (roughly) a $\chi^2$ distribution with 1 degree of freedom testing the null hypothesis that the genotype frequencies derive from $\psi$ following Hardy-Weinberg.


