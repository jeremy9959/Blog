---
type: post
layout: single
title: Mitochondria
excerpt: Scattered info on the mitochondrial genome
tags: mitochondria bioinformatics
---

## Basics

The human mitochondrial genome consists of [16569 base pairs](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) encoding 37 genes:

- 13 protein coding genes (related to the energy production pathways for the cell)
- 22 transfer RNA
- 2 ribosomal RNA

The sequence is extremely compact (there is very little non-coding
DNA) and the [mitochondrial genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) is
slightly different from the standard one. One section of the genome
contains the promoter sites; this region (called the D-loop or the DLP
region) is the most variable; it accounts for 1118 base pairs.



The entire sequence was worked out by [Sanger, et. al.](https://www.ncbi.nlm.nih.gov/pubmed/7219534) and published
in Nature in 1981.  It was one of the first complete genomic sequences to be completed.

## Inheritance and the mitochondrial bottleneck

In humans, all mitochondria are inherited from the mother.  Although sperm do carry mitochondria, they
die out in the fertilized egg.  How this happens is an area of study.
See for example Luo, et. al.  [Sperm Mitochondria in Reproduction: Good or Bad and Where do they Go?](https://www.sciencedirect.com/science/article/pii/S1673852713001641?via%3Dihub).

Because mitochondria are inherited from the mother, and do not
participate in recombination, they are in principle subject to
Muller's Ratchet.  How the mitochondrial genome stays stable despite
this is another area of study.  One aspect of the process is the
*mitochondrial bottleneck.* In this process, at the time of the
production of egg cells, the mitochondrial population of the
progenitor cell is amplified and then distributed randomly among the
daughter cells. This creates a cell population containing mitochondria
populations with varying degrees of genetic diversity.  Selection
among these cells then leads to uniform populations among cells that
survive.  In fact, selection isn't necessary; genetic drift can
account for stabilization of the mitochondrial population.

One interesting consequence of the bottleneck is that the mitochondrial DNA of offspring can vary quite a bit from the mother -- for example, the offspring may be nearly homoplasmic from a very rare variant in the mother's mtDNA.

Various candidate versions of the bottleneck, and other factors that prevent the building up deleterious
mtDNA mutations, are discussed in Mishra and Chan,
[Mitochondrial dynamics and inheritance during cell division, development, and disease](https://www.ncbi.nlm.nih.gov/pubmed/25237825).

In general, *heteroplasmy* refers to the condition in a cell where the mitochondria are genetically diverse; otherwise the cell is *homoplasmic.*
A low level of heteroplasmy is typical in all individuals.

Mitochondria duplicate within a cell and then are distributed to the daughter cells when a cell splits.  Much of this process seems to be random,
and how it is regulated is not clearly understood, although there is coordination between the nuclear genome and the mitochondria.

## Mutations

Heteroplasmy, at a low level, seems to be universal in humans.  In Payne, et. al. [Universal Heteroplasmy of human mitochondrial DNA](https://www.ncbi.nlm.nih.gov/pubmed/23077218),
they found that healthy subjects showed >0.2% heteroplasmy with the level varying between tissue types.  By comparing across relatives, they concluded that
most mutations are probably inherited, and that late-in-life metabolic diseases arising from mitochondrial mutations are the result of clonal expansion of the damaged
variants.  Some mutations are also no doubt somatic.

## Mitochondria in cancer

The mitochondrial genome of cancer cells shows higher mutational load (in the sense of higher
heteroplasmy) than healthy cells.  Whether these mutations are part of the disease mechanism or an
incidental consequence of other features of cancer is not clear.  One hypothesis relating
mitochondria to cancer is the "Warburg effect" which proposes that tumor cells use aerobic pathways
rather than the more efficient oxyidative reactions provided by the mitochondria because of
mitochondrial disfunction.  The 2016 article by Tran, Lee, et. al.
[Targeting Cancer Metabolism - Revisiting the Warburg Effects](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4946416/)
reports that it is "untenable to deny the role of mitochondria in tumurogenesis" and states that cancer cell mitochondria
show:

- reduction in DNA
- lower transcription rate
- accumulation of mutations and deletions.

Table 1 of the article shows, for example, that 70% of colorectal cancers, 60% of ovarian cancers, and 60% of breast cancers
show mutations in the mitochondrial genome,  in the D-loop (the regulatory region), the ribosomal genes, or the coding region.

The paper 2007 [Mitochondrial DNA sequence variation in single cells from leukemia patients](http://www.bloodjournal.org/cgi/doi/10.1182/blood-2006-01-011007)
studies heteroplasmy in leukemia patients and normal controls.  The picture that emerges is complicated; the authors write:

> our analysis for mtDNA variation in single cells from leukemia patients revealed an unexpectedly complex pattern;
> the extent of mtDNA alteration varied greatly among  the patients analyzed and was not generally increased in
> the later stages of leukemia or after medical treatment.

They also found somewhat *fewer* mutations among the leukemia patients than the normal controls, and from this they conclude

> The observed complex pattern in leukemia patients, although on average showing fewer mutations
> than that in controls, suggests that mutations appear at a relatively high rate in leukemic cell
> mtDNA and that homoplasmy (uniformity of sequence within an individual cell) and clonal expansion
> (expression of the fixed mtDNA mutation in a substantial proportion of progeny) can be achieved
> over months rather than requiring decades.

As a very rough cumulative summary of their results, they found that among roughly 100 CD34+ and blast cells from a single leukemia patient,
there were 13-20 haplotypes, each haplotype containing 3-5 cells. Among cells from a healthy patient, they found nearly twice as many
haplotypes, and concluded this difference was significant.  It's worth looking at the table directly because there is a lot of variation in the
numbers.

The 2017 paper by Grandhi, et. al.
[Heteroplasmic shifts in tumor mitochondrial genomes reveal tissue-specific signals of relaxed and positive selection](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5886292/)  attempts to identify selective pressure on mitochondrial genomes in cancer.  In the introduction they assert that
"the presence of selection in somatic cells is unsettled and varies by tissue type."  They use TCGA data to study patterns of mtDNA variants
in tumor and matched normal samples, looking for shifts in heteroplasmy as well as mutations unique to the tumor cells.  They claim that their results show that
mitochondrial genomes generally undergo neutral evolution, and that this allows for "dramatic expansions of frameshift and nonsense" mutations in tumor cells.






