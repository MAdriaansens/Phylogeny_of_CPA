All starting data is publically available data. 
Resources used include HMMs and sequences from Pfam and proteomes/metadata from GTDB. 
NCBI proteomes/genomes was used to generate most of the Eukaryotic database.

On 16th of April 2025 GTDB released its latest version GTDB 226.
On 5th May 2025 we started downloading this data, this included:
         - all metadata from Archaea and Bacteria
         - the proteomes from all representatives

After succesfull download all data will be compiled into a proteome database and tsv for both Archaea and Bacteria
After which this homology searches will be performed

Explaining each directory:

**Database creation
**
**GTDB.**
GTDB proteomic data of representatives was downloaded using -wget from GTDB-Downloads (https://gtdb.ecogenomic.org/downloads).
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_proteins_aa_reps_r220.tar.gz 

_(we recommend doing the downloading and unzipping in a bash script). 
_
We also downloaded the meta data for ar53 and bac120, this contains taxonomic data, completeness etc and will be used to generate the proteomic databases.

A .faa database was created first after which a .tsv file was created matching the protein id with the original id, sequence, gtdb_id, and gtdb taxonomy.
For Archaea and Bacteria the proteome files were seperated using Database_creation/Subset_Fasta_DB.py and Chunk_Bacteria.py. 
For Bacteria the database was subsetted into 61 files (*.faa and *.tsv) while Archaea was subsetted into 12 *.faa files, to allow for parallel processing.

Eukarya data was downloaded by Euk_proteomes_download.py using NCBI ftp links, available in Database_creation/Eukarya/Eukararya_db.tsv.xlsx
We generated, for each taxonomic domain, a fasta file with all the sequences and a .tsv file containing the original protein id, the species id and species name was made.
(done using Database_creation/Eukarya/Create_EukDB.py, Create_TSV_proteinDB_fromGTDB_data.py)


**Pfam.**
Pfam HMMs were downloaded manually and sequences of each pfam were downloaded using a perl script as provided by Pfam.
Code is available on Database_creation/Pfam. After creation of inhouse HMMs these HMMs were concatenated to the PFAM HMMdb and HMMpressed so it can be used for HMMscan in the future.

