**Downloading and generating a Bacteria/Archaea proteomic database from GTDB info**

_GTDB protemes and Metadata downloaded on 15th of May 2025_

**Repeat this step for Bacteria and Archaea**


**Step 0.** 
Download the proteomes of GTDB representatives.
_Do this in a bash script._

    wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_proteins_aa_reps_r220.tar.gz | tar -xzvf gtdb_proteins_aa_reps_r220.tar.gz 
    #followed by gunzip -r only gtdb files.
    #repeat for metadata of Archaea and Bacteria
    
**Step 1.**
On the Proteomic files of the representatives run Create_TSV_and_proteinDB_fromGTDB_data1.py
Input is the proteomic files and the metadata.tsv

_In theory you could stop here but with the size of Bacteria and Archaea databases in mind we recommend splitting the database in pieces._

!! both scripts we use ask for a ceratin size of the chunks make sure they are the same in step2 and step3 !!

**Step2.**
Run Subset_Fasta_DB2.py
Takes as input the *.fasat proteomic database and spits out subsets.
_You can determine the number of subsets and the name of the subsets yourself._
#this is where we stopped for Archaea but we decided to also chunk the .tsv file of BACTERIA



**Step3.**
Run chunk_bacteria.PY
Takes as input the .tsv database and spits out subsets
_You can determine the number of subsets and the name of the subsets yourself._
