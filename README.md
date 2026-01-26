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
Database_creation: 


GTDB data was downloaded using -wget from GTDB-Downloads (https://gtdb.ecogenomic.org/downloads)
Eukarya data was downloaded by Euk_proteomes_download.py using NCBI ftp links, available in Database_creation/Eukarya/Eukararya_db.tsv.xlsx

We generated, for each taxonomic domain, a fasta file with all the sequences and a .tsv file containing the original protein id, the species id and species name was made.

(done using Database_creation/Eukarya/Create_EukDB.py, Create_TSV_proteinDB_fromGTDB_data.py)

For Archaea and Bacteria the proteome files were seperated using Database_creation/Subset_Fasta_DB.py. 
For Bacteria this was subsetted into 61 files (*.faa and *.tsv) while Archaea was subsetted into 12 *.faa files, to allow for parallel processing

Pfam HMMs were downloaded manually and sequences of each pfam were downloaded using a perl script as provided by Pfam.
Code is available on Database_creation/Pfam
