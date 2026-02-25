**Homology search for CPA and IT homologs**

This directory contains 2 subdirectories (for Pfam and Inhouse HMM), as well as scripts for Step1, Step5, Step6 and Step7. 
Step2 is a bash command, Step3 and Step4 are repeats of Step1s MMseqs2 search. 
The search methods differ between CPA and IT. 


**Homology search for IT is 7 steps:**  


**CPA is identical but in Step 1 aside from the HMMsearch using Pfam and the MMseq search using Pfam, we also run a HMMsearches using three CPA specific HMMs. THis the Step1* commandline script**

_Generation of the inhouse HMM is outlined in Step2* HMM generation_
 
  
  **Step 0.** Download Pfam HMM and Pfam sequences (as detailed in Pfam). 

  **Step 1.** Run MMseqs and HMMsearch with the Pfam HMM and HMM sequences as query against all three databases.Retrieve the sequences hit using getting_fasta_from_hit.py
   _HMMsearches must contain the --tblout argument_  
  #the .. script retrieves the protein id and usign that parses the protein.tsv file and writes a fasta file containing the protein id (plus taxonomy) and the sequence.
  #the input for this scripts are:
  1) the MMseq/HMMsearch .tsv file
  2) the protein .tsv file containing the sequences of the sequences you have hit
  3) Wether the .tsv file is MMseq or HMMsearch
  4) the name of your output

#due to eukarya database tsv and bacteria/archaea database tsv being formatted differently a different getting_fasta_from_hit.py is needed when running against Eukarya or Prokarya.

    python getting_fasta_from_hit.py PF03600_MMseqs_e03_vsEukarya_initialsearch.tsv MMSEQ Euk_db_7April_protein.tsv PF03600_MMseqs_e03_vsEukarya_initialsearch.fasta
  **Step 2.**  Concatenate the Bacterial, Eukaryotic and Archaeal sequences into one, using: 

    cat */*/PF03600*fasta > PF03600_vs_all_mergedHMM_mergedMMseq.fasta
  **Step 3** Run another MMseq search using the merged file as query and each database as subject. Retrieve these sequences per database.
  _(this can be done by adding getting_fasta_from_hit.py)_
  

  **Step 4.** Per database run the retrieved against the database again (iterative MMseqs2 search). Retrieve these sequences.

 **Step 5.** Run a quality control by aligning the sequences against their respective HMM, so for PF03600:
          
    hmmalign PF03600_Euksequences_vsPF03600.sthk PF03600 PF03600_Euksequences_iterativeMMseqs_Eukarya.faa
          
  #then we use a script to retrieve the aligned sequence called parse_stockholm.py. The script only retrieves the part of the sequences which has aligned with the HMM (this will be used for HMMscan later on).
  It requires three inputs.
  1) the stockholm file you want to parse
  2) the name you want to give to the output
  3) the 70% cut off, so if a sequence does not align for more than 70% it gets passed. For PF03600 the HMM is 336 aa so the threshold is 236.
     
    python parse_stockholm.py F03600_Euksequences_vsPF03600.sthk F03600_Euksequences_vsPF03600aligned.fasta 236
  **Step 6.** Run HMMscan using the aligned part of the sequences instead of the full length protein. 

  **Step 7.** The sequences whom pass both the HMMscan and HMMalign are the high confidence homologs. 
  Using a python script Step2_homology_search/Parse_All_antiporter_hits_GTDB226Bacteria(H9oB).py.
  We now parse both the HMMscan and HMMalign and retrieve only the sequence ids who pass both. Using the unique taxonomy id attached to the protein id we now also generate a file which contains for each entry in GTDB/In house Eukaryotic database how many CPAs/NhaB/NhaC and NhaD they contain. This script then outputs the PF00999 HMMaligned sequences of all sequences who are high confidence CPA, this will be used later down the line in tree inference.
  
  While I did it differently I added a bit of code so you can retrieve the fl directly from this cript. 
  [In a later iteration I also did the same thing for the fl of the NhaB/NhaC/NhaD]

  It is also possible to run a script called Parse_stockholm_fl_filter.py which does what the original does but retrieves the full length sequences. 
