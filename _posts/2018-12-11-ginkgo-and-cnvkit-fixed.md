---
type: posts
layout: single
excerpt: Notes on [ginkgo](http://github.com/robertaboukhalil/ginkgo.git) and [cnvkit](http://github.com/etal/cnvkit.git)
title: Ginkgo and CNVkit
tags: bioinformatics
---

Algorithms for estimating copy number from low-coverage NGS whole-genome sequencing start from a naive hypothesis that the sequenced fragments
captured during sequencing are selected uniformly at random from the genome.  Therefore, if a region of the genome is repeated n times in a particular
sample, then the number of fragments arising from that region should be n times what one would expect from uniform coverage.

To identify such regions, the typical algorithm:

- subdivides each chromosome in the genome into "bins"
- allocates the reads to each bin
- compares the actual number of reads to the expected number given uniform coverage
- uses a segmentation algorithm to fit a 'square wave' function to the the reads per bin

There are many complications in this process.  For example, it's known that the uniform coverage hypothesis is false:

- Some regions in the genome are not mappable, and so one typically gets no reads there.  Other regions are highly repetitive and the confidence of a particular alignment in that region is low.
- There are biases in the sequencing process, so for example there is a weak but important correlation between the relative abundance of GC vs AT and the number of reads in a given region.

These complications can be addressed, in part, by different ways of constructing the "bins" that are used to identify copy number changes.

Two popular tools for identifying copy number variation in DNA
samples are [ginkgo](http://github.com/robertaboukhalil/ginkgo.git)
and [cnvkit](http://github.com/etal/cnvkit.git).  They take different approaches to the construction of bins.  In subsequent posts we will compare these approaches.

