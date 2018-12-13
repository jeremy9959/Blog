---
type: posts
layout: single
excerpt: How bins are constructed in CNVkit
tags: bioinformatics
title: Binning in CNVkit
---

This discussion pertains to the use of CNVkit to whole genome sequencing (option ```-m wgs``` in the main workflow).

CNVkit takes a different approach to construct bins.  Like ginkgo, it identifies unmappable areas in the genome,
but rather than incorporating those areas into bins, it builds a map of regions around those areas.  More specifically,
it fixes a threshold k and then builds a bed file of regions.  Each region has the property that the longest consecutive
sequence of "N" labels in the genome within that region is shorter than k base pairs.

Typical values of k might be 5000 or 10000 bases.  

Next, having constructed a map of these regions, cnvkit uses one of two methods to subdivide those regions into bins:
- using a targeted average size given as a command line argument.
- estimating a bin size based on a desired number of hits per bin. 

The first of these two methods is done by iterating over the regions in the map of accessible regions.  If the interval is of length
M, and the target size is T, then  if M>=T the region is divided up into bins of size [M/T] or [M/T]+1. Otherwise the interval
is left alone.

The file containing the map of accessible regions is actually quite short. The command to build it (joining over 10000 N's at most) is:

```
$ cnvkit.py access -s 10000 -o access_file.bed reference.fasta
```

Note that ginkgo only cares about the longest stretch of N in each chromosome, arising from the centromere.

