---
type: posts
layout: single
tags: bioinformatics
title: Some fundamental terminology on  time series
excerpt: Stationarity, autocorrelation, ergodicity, ...
---

A time series is a sequence $X_n$ of random variables indexed by the integers $\mathbf{Z}$.  It is a particular case
of a stochastic process.  

The time series is *strictly stationary* if, for any integers $m\ge 0$ and $k>0$ we have that
the joint distributions of $(X_0,\ldots, X_m)$ and $(X_k,\ldots, X_{m+k})$ are the same. (This is the definition
given in Resnick.  Wikipedia gives a stronger condition, in which one requires the joint distribution of 
$X_t$ at a finite, but arbitrary set of times, to be invariant under shifts in time. )

Notice that a time series that is zero when $n<0$ and nonzero when $n\ge 0$ isn't strictly stationary. 

The time series is *wide-sense stationary* if the function $m(t)=E(X_t)$ is time-invariant, the
second moments $E(|X_t|^2)$ are finite for all $t$, and the correlation between $X_t$ and $X_{t'}$ depends
only on the time interval between $t$ and $t'$.

A wide-sense stationary process  is "mean-ergodic" if the random variable 
$$
\hat{\mu}(N)=\frac{1}{N}\sum_{1}^{N} X_t 
$$
converges to the (constant) mean $E(X_t)$ in $L^{2}$ as $N\to\infty$.  In other words,  $E((\hat{\mu}(N)-\mu)^2)\to 0$ as
$N\to\infty$.  When this holds, one may estimate the process mean from time averages.

From these [course notes](http://ece-research.unm.edu/bsanthan/ece541/ergmean.pdf) we see the relationship between
ergodicity in the mean and the autocovariance.  We have

$$
Var(\hat{\mu}(N))=E((\hat{\mu}(N)-\mu)^2)
$$

$$
=E\left\{\left(\frac{1}{N}\sum_{1}^{N} X_{t_1} - \mu\right)\left(\frac{1}{N}\sum X_{t_2}-\mu\right)\right\}
$$

Using the definition of the (auto)-covariance, this reduces to

$$
Var(\hat{\mu}(N))=\frac{1}{N^2}\sum_{1}^{N}\sum_{1}^{N} E((X_{t_1}-\mu)(X_{t_2}-\mu))
$$

The expectation term inside the sum is the definition of the autocovariance $C(t_1,t_2)$ which by the stationarity
hypothesis only depends on $|t_1-t_2|$.  

Simplifying the calculation by using integrals instead of sums, the notes quoted above ultimately conclude
that the $L^2$ condition amounts to

$$
\lim_{T\to\infty}\frac{1}{T}\int_{-T}^{T}C(z)(1-\frac{|z|}{T}) dz =0.
$$

Note that if the autocovariance function is compactly supported this is automatically satisfied.  Otherwise it's
a growth condition on that function -- it's sufficient if the area under $C(z)$ above $[-T,T]$ grows slower than $T$ --
but in fact the condition is a bit weaker than that because the weight function is a triangular bump of height 1 and base
$[-T,T]$. 

