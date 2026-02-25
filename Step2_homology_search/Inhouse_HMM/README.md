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
It is important to note that we cannot use HMMscan and HMMalign during curation of sequences, since we are working under the assumption that the HMMs are not perfect for detecting distant sequences. So early usage of HMMs would remove these distant sequences, therefore an e-value of 1*10^-5 was chosen. 

We then ran hmmsearch using each hmm, and the retrieved sequences were then aligned against the Pfam00999, using HMmalign, and checked how many of those mapped over 70% to the hmm.

The HMM with the highest number of sequences was then chosen per taxonomic domain. If one domain had multiple HMMs retrieving the same amount we chose the HMM generated using the most sequences as this would be a HMM representing more diversity.


**Comparing HMMs**
Finally we compared the number of sequences the 3 HMMS retrieved after running a HMmalign 70% filter.

_It was shown that the PF00999 HMM was often outperformed when ran against each taxonomic domain compared to each HMM,  but never the worst HMM._


_Parameters selection (detailed) and steps_
**Step 1. a MMseqs2 was run,** it was observed from early attempts that a HMM always performed better when it was derived from a MMseq search with a coverage threshold was set in place. In the beginning we also tested e-value but observed that e-value 1*10^-5 performed best, since we do not have a QC filter inplace (such as HMMalign or HMMscan). Sequences were then retrieved using getting_fasta_from_hit_extra_Euk.py

    declare -a array=(0.5 0.3)
    Options for MMseq --cov 0.5 / --cov 0.3
    
    module load MMseqs2/15-6f452-gompi-2023a
    mmseqs easy-search -e 1.00E-03 --cov ${array[$SLURM_ARRAY_TASK_ID]} --threads 12 SG_seed_update_Jan14_2024_135equences.faa ${DB} Manual_sequences_vsEukarya_e03_cov${array[$SLURM_ARRAY_TASK_ID]}.tsv tmp
    module purge
    module load Python/3.11.6-foss-2023a
    python getting_fasta_from_hit_extra_Euk.py Manual_sequences_vsEukarya_e03_cov${array[$SLURM_ARRAY_TASK_ID]}.tsv  MMSEQ ${Euk_TSV} Manual_sequences_vsEukarya_e03_cov${array[$SLURM_ARRAY_TASK_ID]}.fasta
**Step 2. Cluster retrieved sequences,** Now we got a large set of sequences using our manual selected set and MMseqs2. Running ** Sequence_clustering.sh** we could run an array of MMseqs2 commands on one or mutliple sequences. We run this clustering step to reduce redundancy but keep the more distant sequences. 

    sbatch Sequence_clustering.sh
