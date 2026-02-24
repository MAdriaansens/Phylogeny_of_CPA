In this directory are the CPA specific HMM per Taxonomic domain.
They were chosen based on the number of sequences used in their alignment (for best representation)
As well as the number of sequences passed the HMMalign when a HMMsearch was performed with these HMMs.

!!The Eukarya HMM is 'Manual_seq_cov30_e050_seqid0.7_genafpair_aligned.hmm!!
**
**The premise of generating an inhouse HMM is that for CPA the diversity is still largely unmapped.**

_Meaning that any HMM generated using characterized sequences will be biased._



To resolve this we generated 3 HMMs per taxonomic domain. Using manually selected seed sequences (SG_seed_update_Jan14_2024_135equences.faa) as the starting query.

These sequences were by Pfam and UniProt annotated as CPA or were identified by blasting Halophilic proteomes against CPA representatives.



_The generation of each HMM is detailed in 'Supplementary methods'._

**HMM selection and parameter selection**

Since we do not exactly know which parameters for HMM generation were best we used and tested out most. 

We then ran hmmsearch using each hmm, and the retrieved sequences were then aligned against the Pfam00999, using HMmalign, and checked how many of those mapped over 70% to the hmm.

The HMM with the highest number of sequences was then chosen per taxonomic domain. If one domain had multiple HMMs retrieving the same amount we chose the HMM generated using the most sequences as this would be a HMM representing more diversity.


**Comparing HMMs**
Finally we compared the number of sequences the 3 HMMS retrieved after running a HMmalign 70% filter.

_It was shown that the PF00999 HMM was often outperformed but never the worst HMM._
