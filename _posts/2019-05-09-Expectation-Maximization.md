---
layout: single
type: post
title: Expectation Maximization
excerpt: Proof of the EM algorithm 
tags: math bioinformatics
---
## Expectation Maximization

### Introduction
I've had a lot of trouble finding a comprehensible explanation of expectation maximization, even in simple cases like those that arise when trying to reconstruct haplotype frequencies from phenotypes.  So here is my attempt.

The relevant situation is when you have some observed data which derives from hidden data and you have a statistical model, depending on some parameters, that predicts the hidden data and hence the observed.  You'd like
to estimate the parameters of the model using only the observed data.

This situation comes up, for example, with hidden Markov models, but an even simpler situation is the fairly typical problem from genetics that is
discussed in [this post](https://jeremy9959.github.io/Excoffier-Slatkin).  That post considers the situation
in which you have two loci, each of which has two alleles, so that there are four haplotypes possible.
Looking at an individual, you can observe each locus separately but you can't directly observe the haplotypes.
So for example, an individual with phenotype AaBb could have haplotype AB/ab or Ab/aB.  If one assumes a model where the two loci are independent and mating is random -- so that phenotypes are found by choosing a pair of haplotypes from the pool at random -- one can use EM to estimate the most likely haplotype frequencies.

Using likelihood ratios one could go further and test the model to see how well it fits the observed data.

### Setup

We assume that we have a random variable $X$ depending on a statistical model with parameters $\theta$.
We assume further that there is an additional random variable $Z$, also depending on $\theta$, that we cannot observe
directly, but which in some way is connected to the values of $X$.  We'd like to find the maximum likelihood estimate for $\theta$ given observations of $X$.

Having observed $x_1,\ldots, x_N$, we consider the log-likelihood of $\theta$ 
$$
\ell(\theta)=\log P(\{x_{i}\}|\theta)=\sum_{i=1}^{N}\log P(x_{i}|\theta)
$$
where we've implicitly assumed an uninformative prior on $\theta$ and that the $x_i$ are chosen $iid$.

We assume that the $x_i$ are derived from the hidden $Z$ so
$$
P(x_i|\theta)=\sum_{\alpha} P(x_i,z=\alpha|\theta)
$$
where we've recovered $P(x|\theta)$ by integrating out the possible values of $z$.
So purely formally we have
$$
\ell(\theta)=\sum_{i=1}^{N}\log\sum_{\alpha}P(x_i,z=\alpha|\theta).
$$

#### Jensen's Inequality

Theorem: 
Let $\phi$ be a convex function on the real line and let $X$ be an integrable real-valued random variable on a probability space $\Omega$.  Then $E[\phi(X)]\ge \phi(E[X])$.  If $\phi$ is strictly convex, then equality holds only if $X$ is ae constant.

Notice that $-\log(x)$ is strictly convex, so in particular
$$
-\frac{1}{m}\sum_{i=1}^{m}\log(x_i)\ge -\log(\frac{1}{m}\sum_{i=1}^{m} x_{i})
$$
which says that the geometric mean is less than the arithmetic mean, so in some sense Jensen's inequality is a vast generalization of that fact.




### Applying Jensen's Inequality

  We take advantage of Jensen's inequality as follows.  Suppose we put a probability distribution on the $z$,  so that
we have, for each value $\alpha$ of $z$, and each $i$, a probability $Q^{(i)}_{\alpha}=Q(z=\alpha)$ with $\sum_{\alpha} Q^{(i)}_{\alpha}=1$.  Trivially, we have
$$
\ell(\theta) = \sum_{i=1}^{N}\log\sum_{\alpha}Q^{(i)}_{\alpha}\frac{P(x_{i},z=\alpha|\theta)}{Q^{(i)}_{\alpha}}
$$
Now viewing $P(x_i,z=\alpha|\theta_0)$ as a function of $\alpha$ (for a given $x_i$), the inner sum is 
$$
\sum_{\alpha}Q^{(i)}_{\alpha}\frac{P(x_i,z=\alpha|\theta)}{Q^{(i)}_{\alpha}}=E_{Q^{(i)}}[\frac{P(x_i,z=\alpha|\theta)}{Q^{(i)}_\alpha}]
$$
where the expectation is taken over the probability distribution $Q$. 

Using Jensen's inequality, we have
$$
\log E_{Q^{(i)}}[\frac{P(x_i,z=\alpha|\theta)}{Q^{(i)}_\alpha}]\ge \sum_{\alpha} Q^{(i)}_{\alpha}\log\frac{P(x_i,z=\alpha|\theta)}{Q^{(i)}_\alpha}
$$

Now we have
$$
\ell(\theta)\ge \sum_{i=1}^{N} \sum_{\alpha} Q^{(i)}_{\alpha}\log\frac{P(x_i,z=\alpha|\theta)}{Q^{(i)}_\alpha}=\sum_{i=1}^{N} \sum_{\alpha} Q^{(i)}_\alpha\log P(x_i,z=\alpha|\theta) - \sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha\log Q^{(i)}_\alpha
$$

### The optimization problem

Now let's imagine that we have an initial guess for $\theta$.  Call this guess $\theta_0$.  Our goal
is to improve our estimate of the maximum likelihood by finding a $\theta$ so that $\ell(\theta)>\ell(\theta_0)$.  Define
$$
Q_{\alpha}^{(i)}(\theta)=P(z=\alpha|x_i,\theta)=\frac{P(x_i,z=\alpha|\theta)}{\sum_{\alpha}P(x_i,z=\alpha|\theta)}=\frac{P(x_i,z=\alpha|\theta)}{P(x_i|\theta)}
$$

Observe first that, tracing back through Jensen's inequality, we have
$$
\ell(\theta_0)= \sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha(\theta_0)\log P(x_i,z=\alpha|\theta_0) - \sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha(\theta_0)\log Q^{(i)}_\alpha(\theta_0) = \sum_{i=1}^{N}\log P(x_i|\theta_0)
$$

Now we look at the terms on the right side of the inequality for $\ell(\theta)$.  We view the probabilities inside the logarithms as functions of $\theta$, but we set the outside $Q_\alpha^{(i)}(\theta)$ to their values for $\theta=\theta_0$.  Let's first look at the term
$$
H(\theta,\theta_0)=\sum_{i}\sum_{\alpha}Q_{\alpha}^{(i)}(\theta_0)\log Q_{\alpha}^{(i)}(\theta).
$$
For any fixed value of $\theta$ and $i$, both $Q_{\alpha}^{(i)}(\theta)$ and $Q_{\alpha}^{(i)}(\theta_0)$ are probability distributions over $\alpha$ and hence sum to $1$.  [Gibb's inequality](https://en.wikipedia.org/wiki/Gibbs%27_inequality) says that if $q_{\alpha}$ and $p_{\alpha}$ are any two probability distributions over $\alpha$, then 
$$
\sum_{\alpha} q_{\alpha}\log p_{\alpha}
$$
takes its maximum value exactly when $p=q$. Consequently, if
$$
\sum_{i=1}^{N} \sum_{\alpha} Q^{(i)}_\alpha(\theta_0)\log P(x_i,z=\alpha|\theta) >
\sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha(\theta_0)\log P(x_i,z=\alpha|\theta_0)
$$

then

$$
\sum_{i=1}^{N} \sum_{\alpha} Q^{(i)}_\alpha(\theta_0)\log P(x_i,z=\alpha|\theta)- \sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha(\theta_0)\log Q^{(i)}_\alpha(\theta)
$$

is greater than

$$
 \sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha(\theta_0)\log P(x_i,z=\alpha|\theta_0) - \sum_{i=1}^{N}\sum_{\alpha}Q^{(i)}_\alpha(\theta_0)\log Q^{(i)}_\alpha(\theta_0)=\ell(\theta_0)
$$

The conclusion of this computation is that, to find $\theta$ with $\ell(\theta)>\ell(\theta_0)$, it suffices to find the value of $\theta$ that maximizes the function

$$
\Theta(\theta,\theta_0)=\sum_{i}\sum_{\alpha} Q_{\alpha}^{(i)}(\theta_0)\log P(x_i,z=\alpha|\theta)
=\sum_{i}\sum_{\alpha} P(z=\alpha|x_i,\theta_0)\log P(x_i,z=\alpha|\theta)
$$

### Iteration

After finding the maximum $\theta$, we repeat the entire computation with $\theta_0$ set to this new $\theta$.

### More on the maximization step

What makes the EM algorithm attractive is that it replaces the complicated likelihood with the potentially simpler problem of maximizing
$$
\Theta(\theta,\theta_0)=\sum_{i}\sum_{\alpha} P(z=\alpha|x_i,\theta_0)\log P(x_i,z=\alpha|\theta).
$$

For example, if $P(x_i,z=\alpha)$ is a normal distribution then its logarithm is linear.  There are a number of places where such an example is worked out.

In simple genetics situations it's common for $P(x_i,z=\alpha)$ to be multinomial.  In the prototypical case, suppose we are looking at Mendelian genetics and we have a dominant and recessive allele with population frequencies $p$ and $q$ with $p+q=1$.  

Let's make an initial guess $p_0$ and $q_0$ for the frequencies. 

We can observe the number $D$ of phenotypes containing at least $1$ dominant allele, and $R$ the number of homozygous recessive individuals.  There are three genotypes (these are the values of $z$).   We have $P(x=D,z=AA)=p^2$, $P(x=D,z=Aa)=2pq$, and $P(x=R,z=aa)=1$.  The other possibilities for $P(x,z)$ are zero.  

Now 
$$P(z=AA|x=D,p_0,q_0)=p_0^2/(p_0^2+2p_0q_0)$$

$$P(z=Aa|x=D,p_0,q_0)=2p_0q_0/(p_0^2+2p_0q_0)$$

and
$$P(z=aa|x=R)=1.$$

We observe $N_D$ dominant phenotypes and $N_R$ recessive ones.  Then our $\Theta$ function is

$$
\Theta=\frac{p_0^2N_D}{p_0^2+2p_0q_0}\log(p^2)+\frac{2p_0q_0N_D}{p_0^2+2p_0q_0}\log(2pq)+N_R\log(q^2).
$$

or

$$
\Theta = (\frac{2p_0^2N_D+2p_0q_0N_D}{p_0^2+2p_0q_0})\log(p)+(\frac{2p_0q_0N_D}{p_0^2+2p_0q_0}+2N_R)\log(q)+C
$$

where $C$ is a constant.

Since $p+q=1$ it's easy to see that the maximum likelihood estimates for $p$ and $q$ are 

$$
p=\frac{1}{1+q_0}\frac{N_D}{N_D+N_R}
$$

and

$$
q=\frac{(q_0/(1+q_0))N_D+N_R}{N_D+N_R}
$$

so iterating the formula

$$
q'=\frac{(q/(1+q))N_D+N_R}{N_D+N_R}
$$

leads to the maximum likelihood estimate for $q$.  Of course, looking at fixed points we see that
this is just the formula

$$
q_{\mathrm{max}}=\sqrt{R/(R+D)}
$$
