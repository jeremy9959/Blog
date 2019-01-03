---
type: post
layout: single
title: "FPKM and TPM Units for RNA Seq after Pachter"
excerpt: Pachter explains units of expression
tags: rna-seq bioinformatics
---
# FPKM and TPM units for RNA seq

See [Pachter, L. Models for transcript quantification from RNA-seq](https://arxiv.org/pdf/1104.3889.pdf).

## Count-based models

### RPKM/FPKM

Simplest version: transcriptome is a set of transcripts with different abundances, and a read is produced by choosing a site in a transcript for the read uniformly at random from among all positions.

* $T$ is the set of transcripts.
* $\ell(t)$ is the length of transcript $t$ for $t\in T$.
* $\rho(t)$ is the relative abundance of transcript $t$ in the transcriptome.  (That is, it is the number of mRNA molecules of type $t$, over the total number of such molecules).
* $m$ is the length of the reads
* $F$ is the set of all possible reads and $F_t$ is the set of reads that map to transcript $t$. 
* The effective length $\tilde{\ell}(t)$ is $\ell(t)-m+1$, the
total number of places that a read of length $m$ could begin.

What's the chance of choosing a read from a transcript $t$? Let $m_t$ be the number of transcripts of type $t$ and $M=\sum m_{t}$ be the total number of transcripts.  The number of reads coming from 
$m_t$ is $m_t \tilde{\ell}(t)$.  The total number of reads from all transcripts is $\sum_{i\in T} m_{i}\tilde{l}(i).$  Therefore the
probability of getting a read from $t$ is
$$
\alpha_t = \frac{m_{t}\tilde{\ell}(t)}{\sum_{t} m_{t}\tilde{\ell}(t)}
$$
and, since $\rho_t=m_{t}/M$, this is the same as
$$
\alpha_t=\frac{\rho_t\tilde{\ell}(t)}{\sum_{t} \rho_{t}\tilde{\ell}(t)}
$$

It's worth observing that $\sum \alpha_{t}=1$ and the $\alpha$'s are non-negative.

The $\alpha$ and $\rho$ distributions are related (via the lengths $\tilde{\ell}$) in the opposite direction by the equation
$$
\rho_{t} = \frac{(\alpha_t/\tilde{\ell(t)})}{\sum_{i\in T} \alpha_{i}/\tilde{\ell}(i)}
$$

The maximum likelihood estimate for the $\alpha_t$ are $X_t/N$
where $N$ is the total number of mapped reads and $X_t$ is the number of reads mapped to $t$. 

The maximum likelihood estimate for $\rho_t$, which is the relative abundance of transcript $t$ among the expressed transcripts, is
(using the above equations) PROPORTIONAL TO
$$
\hat{\rho}_{t}\sim\frac{X_{t}}{N\tilde{\ell}(t)}.
$$
(You need to divide by the sum of the terms on the right to get equality).

The number $$\frac{X_{t}}{\tilde{\ell}(t)N}(10^{9})$$ is the "Reads per Kilobase per millions of mapped reads" because $$X_{t}/(N/10^{6})$$ is the "Reads per million mapped reads" and then you divide that by $$\tilde{\ell}(t)/10^3$$ to get "Reads per million mapped reads per kilobase."

### TPM

The 'transcripts per million' measure (maximum likelihood) is $10^6\rho_{t}$ so to compute TPM you take the abundance measured in  RPKM and divide by the sum of RPKM over all the transcripts.
