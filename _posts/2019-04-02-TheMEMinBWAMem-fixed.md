---
type: post
layout: single
title: The mem in bwa and what it means
excerpt: MEM means Maximal Exact Match -- and it can trip you up.
tags: bioinformatics
---

## BWA, seeds, and maximal exact match.

The BWA aligner is usually run with the mem option for sequences above
70bp.  As the [manual page](http://bio-bwa.sourceforge.net/bwa.shtm)
shows, there are many options to this program, enabling the user to
tailor the alignment algorithm to specific situations.  One aspect of
the alignment algorithm tripped me up in a particular (artificial)
case, and here is an explanation of that case and what I believe was
happening.

The mem option stands for Maximal Exact Match.  When bwa mem begins
its alignment process for a particular read, the first thing it does
is look for a long substring that matches the reference exactly and
which can't be extended to a longer match at either end -- that is, it
finds a maximal exact match.  If that maximal exact match is too
short, then bwa discards that read (marks it as unmapped) and moves
on.  The threshold for this is called the 'minimum seed length' and it
is set by default at 19.

Suppose, then, that we have, say, a 100 base read with 10 SNPs space
across the read so that the longest match against the reference is 10bp
long.  Then that read won't get mapped to the reference.

That is, unless one uses the -k option to change this minimum seed
length to a smaller number.

I discovered this in a very artificial situation: I generated a "fake"
fastq file by:

- starting with the human mitochondrial genome, 16569 bp long
- picking 1000 random locations and introducing a SNP at those
locations, yielding an alternative genome.
- chopping up that alternative randomly into 1000 76bp strings and choosing
qualities at random

Running ```bwa mem``` on this fastq file with the original
mitochondrial genome as reference was yielding a bunch of unmapped
reads, which seemd strange under the circumstances.

After looking closely at the ```bwa mem``` algorithm I finally
realized that the problem was that my random choice 1000 snps on a
16kb genome yielded roughly one snp per 16 or 17bp; so it wasn't that
unlikely that I'd have a non-trivial number of 76bp reads with a
maximal exact match below the 19 threshold.

Lowering the threshold eliminated the unmapped reads (although
lowering it too much meant that the short reads at the end of the
genome got mapped fairly randomly).


## Re-seeding

There's a second option related to seeds in the bwa man page.  The -r
option, which takes a floating point number FLOAT, is described as follows:
```
Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is
a key heuristic parameter for tuning the performance. Larger value
yields fewer seeds, which leads to faster alignment speed but lower
accuracy.
```
The default value is 1.5.  What does this mean?  A problem for another
day.








