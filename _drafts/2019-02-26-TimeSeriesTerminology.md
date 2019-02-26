---
type: posts
layout: single
tags: bioinformatics
title: Some fundamental terminology on  time series
excerpt: Stationarity, autocorrelation, ergodicity, ...
---

A time series is a sequence $X_n$ of random variables indexed by the integers $\mathbf{Z}$.  It is a particular case
of a Stochastic process.  

The time series is *strictly stationary* if, for any integers $m\ge 0$ and $k>0$ we have that
the joint distributions of $(X_0,\ldots, X_m)$ and $(X_k,\ldots, X_{m+k})$ are the same. (This is the definition
given in Resnick.  Wikipedia gives a stronger condition, in which one requires the joint distribution of 
$X_t$ at a finite, but arbitrary set of times, to be invariant under shifts in time. )

Notice that a time series that, for example, is zero when $n<0$ and nonzero when $n\ge 0$ isn't strictly stationary. 
