---
type: post
layout: single
title: deFinetti's Theorem
excerpt: some notes on deFinetti's Theorem and related topics
tags: math
---

## De Finetti's Theorem

De Finetti's theorem is a fundamental result in Bayesian probability and is closely related to the theory of the
Dirichlet Distribution and the Dirichlet Process which arise in clustering. 

For the first part of this post we follow the lovely paper [An elementary proof of de Finetti's Theorem](https://arxiv.org/pdf/1809.00882.pdf) by [Werner Kirsch](https://www.fernuni-hagen.de/stochastik/team/werner.kirsch.shtml).

**Theorem:** (de Finetti) Let $S=(X_1,X_2,\ldots)$  be an infinite sequence of $0,1$-valued random variables
that form an "exchangeable sequence,", meaning that $P(X_1=a_1,X_2=a_2,\ldots, X_N=a_N)$ is invariant under permutation
of the $X_i$ for any $N$.  Then there exists a probability measure $\mu$ on $[0,1]$ such that, for any $N$ and any sequence
of zeros and ones $a_1,\ldots, a_N$, we have

$$
P(X_1=a_1,\ldots, X_N=a_N)=\int y^e(1-y)^s d\mu(y)
$$

where $e$ is the number of $a$'s equal to $1$ and $s$ is the number of $a$'s equal to zero (so $e+s=N$).

This result is usually interpreted by saying that the infinite exchangeable sequence is a mixture of iid Bernoulli random
variables with mixture measure $\mu$ -- so that to sample from the distribution 
one first "picks" a probability $p$ from the distribution $\mu$ and then does $N$ 
flips of a bernoulli coin with heads probability $p$.

The fact that we have an *infinite* sequence of random variables is crucial for this result.  See for example
Persi Diaconis's early, and elementary, paper 
[Finite forms of de Finetti's theorem on exchangeability](http://statweb.stanford.edu/~sabatti/Stat370/synthese.pdf) for the beginning
of a long story on the relationship between finite and infinite sets of exchangeable random variables.

To get a sense of how this all works, notice first of all that the "exchangeability condition" means that all of the
information in the sequence $S$ can be encoded in numbers $P(k,m)$ for $m\in\mathbb{N}$ and $0\le k\le m$ where

$P(k,m)$ is the probability of getting $k$ ones from a choice of $m$ variables $X_i$ from $S$.


In fact, suppose we have a collection $(a_1,\ldots, a_m)$ of zeros and ones, and $k=\sum a_i$ is the number of ones.  Then

$$
P(X_{i_1}=a_1,X_{i_2}=a_2,\ldots,X_{i_m}=a_m) = \frac{P(k,m)}{\binom{m}{k}}.
$$

since exchangeability means that the location of the $k$ $1$'s among the $m$ slots is uniformly distributed among the
$\binom{m}{k}$ possible sites.

**Lemma:** The probability measure $\mu$ defined by the function $P(k,m)\to \mathbb{R}$ described above
is uniquely determined by the values $P(n,n)$.

**Proof:** The fact that the $P(k,m)$ determine a measure on the full product space means that

$$
\begin{array}
P(X_1=a_1,\ldots,X_n=a_n)& = & P(X_1=a_1,\ldots,X_n=a_n,X_{n+1}=0) \cr
& & +P(X_1=a_1,\ldots,X_n=a_n,X_{n+1}=1).
\end{array}
$$

Therefore

$$
\frac{1}{\binom{n}{k}}P(k,n)=\frac{1}{\binom{n+1}{k}}P(k,n+1) + \frac{1}{\binom{n+1}{k+1}}P(k+1,n+1).
$$

Rearranging and simplifying we obtain the relation

$$
P(k,n+1) = \frac{n+1}{n+1-k}P(k,n) - \frac{k+1}{n+1-k}P(k+1,n+1).
$$

From this relation, if we know $P(k,n)$ for all $n<N$ and $P(N,N)$, we can recursively compute
$P(k,N)$ for $0\le k\le N-1$ to obtain all $P(k,N)$.

Notice also that the random variable

$$
S_N = \frac{1}{N}\sum_{i=1}^{N} X_i
$$

counts the proportion of $1$'s among the first $N$ of the $X_i$, so its distribution is governed by the $P(k,m)$,
so that for $0\le k\le N$ we have $S_N$ supported on the fractions $k/N$ and

$$
P(S_N=k/N) = P(k,N).
$$

## Kirsch's Proof

The strategy of the proof is to exploit some theoretical results relating convergence of moments to weak convergence of measures.
This is a somewhat subtle point, since the polynomials don't belong to
the bounded continuous functions so can't be used directly as "test functions" for weak convergence.  

**1.** If $\nu$ is any probability measure on $\mathbb{R}$, and $\mu$ is a probability measure with support contained in $[a,b]$,
and if $\nu$ and $\mu$ have the same moments, then $\nu=\mu$.  

**Proof:** This is [Theorem 2.55](https://www.fernuni-hagen.de/stochastik/downloads/momente.pdf) on page 28 of Kirsch's book
on moments.  That theorem asserts that if $\mu$ is a bounded measure with "moderately growing moments" and $\nu$ is a bounded
measure $\nu$ such that $\mu$ and $\nu$ have the same moments, then they are equal.  Moderate growth means that the $k^{th}$
moment is bounded by $AC^{k}/k!$ for constants $A$ and $C$ and all $k$.  If $\mu$ is compactly supported, as in our case,
then this condition is automatic.  The proof compares the characteristic functions (Fourier transforms)
of $\mu$ and $\nu$.  The proof also relies on Prohorov's theorem.  

**2.** If $\mu_n$ form a sequence of probability measures on $[0,1]$
all of whose moments $m_k(\mu_n)$ converge to some $m_k$ as $n\to\infty$, then the $\mu_n$ converge weakly to a unique probability
measure on $[0,1]$ with moments $m_k$. 

**Proof:** This is a consequence of [Theorem 2.56](https://www.fernuni-hagen.de/stochastik/downloads/momente.pdf) on page 28 of Kirsch's book (which is quite a bit more general).  The sequence of measures $\mu_n$ has bounded moments (even bounded second moment is enough),
so it is *tight*.  Since all are probability measures, we know that the sequence $\mu_n$ has a weakly convergent subsequence by
Prohorov's theorem. Take any such convergent subsequence and let $\mu$ be the limit of that subsequence.  The moments of $\mu$ are the limits $m_k$.  Any other convergent subsequence has a limit $\mu'$ which has the same moments $m_k$, so $\mu=\mu'$.  In other words,
every convergent subsequence of the $\mu_n$ has the same limit, so the sequence $\mu_n$ converges (weakly) to $\mu$.

**Proposition:** The random variables $S_N$ converge in distribution to a probability measure $\mu$ on $[0,1]$ as $N\to\infty$
with moments
$$
\int y^{k}d\mu(y) = P(k,k).
$$

**Proof:** Given **1** and **2** above, we need to check that the moments of the distribution $\mu_N$ determined by $S_N$ converge
to $P(k,k)$ as claimed.  In other words, we need to compute
$$
\lim_{N\to\infty} E((\frac{1}{N}\sum_{i=1}^{N} X_i)^{k}) = \lim_{N\to\infty}\frac{1}{N^{k}}\sum_{i_1,\ldots,i_k} E(X_{i_1}X_{i_2}\cdots X_{i_k})
$$
Since the individual random variables $X_i$ are idempotent, and because of exchangeability, to work this out  we need to count
how many times each monomial equivalent to $X_1 X_2\cdots X_r$ for $1\le r\le k$ occurs in the sum. For case where all
$X_i$ are distinct, we have the number of ways to choose a subset of size $k$ from $N$ elements, or $\binom{N}{k}$; and then
each subset contributes $k!$ terms.  Thus the number of terms of this form grows with leading term $N^{k}$.  

For all of the other situations, we need to count subsets of size $k$ chosen from $N$ elements with *at least one repetition.*
To estimate the number of these, we choose at most $k-1$ indices between $1$ and $N$ for a total of choices on the order
of $N^{k-1}$; and then we make a list of length $k$ of choices from this set with $k-1$ elements, for a total of at most
$(k-1)^kN^{k-1}$ choices.  (In fact this is a massive overcount but we only care about the leading order of $N$).  

It follows that, as $N\to\infty$, only the terms with $k$ distinct $X_i$ survive; and in the limit we obtain

$$
\lim_{N\to\infty} E((\frac{1}{N}\sum_{i=1}^{N} X_i)^{k}) = E(X_1 X_2\cdots X_k).
$$

Finally, the product $X_1 X_2\cdots X_k$ is zero unless all of the $X_i$ are $1$, and thus this expectation is exactly $P(k,k)$.

This proof tells us that the De Finetti measure $\mu$ is characterized by the $P(k,k)$:

$$
\int y^{k} d\mu = P(k,k).
$$

To complete the proof, 
observe that if we set

$$
 \frac{1}{\binom{m}{k}}P(k,m) = \int y^{k}(1-y)^{m-k} d\mu ,
$$

where the binomial coefficient appears since we need to consider all possible ways of getting $k$ ones in $m$ tries,
the identity

$$
\int y^{k}(1-y)^{m+1-k}d\mu = \int y^{k}(1-y)^{m-k} d\mu - \int y^{k+1}(1-y)^{m-k}d\mu
$$

so probabilities determined by this limiting $\mu$ satisfy the same recurrence as the $P(k,m)$ and agree on $P(k,k)$
and are therefore the same.








