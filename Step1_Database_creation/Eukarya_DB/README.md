**Download and generate the Eukaryotic proteomic database**

_All proteomic files were downloaded on 16th of May 2025._

The input for the creation of the eukaryotic database is the 'Eukarya_metadata.tsv' file.
This file contains the species information, taxonomy and more of 245 species including ftp links to NCBI.
This dataset is custom generated and is informed by work of Eva van Deutekom (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007301#sec008) 
As well as Carlos Santana-Molina (https://www.nature.com/articles/s41559-025-02648-0).

**Step 1**

  The Euk_proteomes_download1.py script basically reads the Eukarya_db.tsv.xlsx.
    It takes the ftp link and runs 'wget' to download sequences. It also spits out a list of species who need JGI downloads. 
    
For some species NCBI does not contain proteomic data but JGI or other data sources might (i.e. the JGI download list). Thus do not expect the downloading to be completely automated.
    We strongly suggested using wc -l to assess presence of all 245 species.

    ls -lh */* | wc -l

After you have made sure you have downloaded all 245 files (be aware some might be zip files and or in subdirectories) we start step 2.

**Step 2**

Actually creating and naming the proteomic files. 
This is done using Create_EukDB2.py, which takes the downloaded proteomes.
It also generates two files: 
  1. A protein database of all eukaryotic sequences (with sequence ids agnostic of species (Euk))
  2. A tsv file containing taxonomic data, species id, protein id and sequence (feel free to edit the name of the protein how ever you want)

**After running both scripts and doing some testing you should have:**

245 proteomes.
A database of eukaryotic proteomes with sequence ids agnostic of species origin
A .tsv file containing the new protein id, taxonomy, and sequence (among other information)
