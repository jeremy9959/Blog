---
type: posts
layout: single
tags: bioinformatics
title: Amplification Bias - 2
excerpt: See Calibrating genomic and allelic coverage... by Zhang, et. al.
---

This is a further look at the paper:
[Calibrating genomic and allelic coverage bias in single-cell sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4922254/) by Zhang, et. al.

In particular we reproduce Figure 2a cited in our [earlier post](/AutoCorrelation).  The data set can be retrieved
as a bam file from SRA reference [SRR1779331](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1779331).  The coverage
data can be extracted from chromosome 1 (which is the only contig that's used for the graph) using ```samtools depth.```

To compute the autocorrelation, after considering several techniques, we settled on exploiting the sparse matrix library
from scipy.  That is described in [this post](/ComputingAutoCorrelation).

The figure we obtained does in fact reproduce Zhang:

<img src="{{ "/assets/images/zhang2a_repro.png" | relative_url }}">


