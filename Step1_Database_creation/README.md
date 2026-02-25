This is the code used to create the Archaea DB, Bacteria DB and Eukarya DB.

Files used but missing are the metadata files from GTDB from Archaea and Bacteria respectively as well as the Pfam HMM DB.
These files are available publically. 

Run the two json scripts for Eukarya, Archaea and Bacteria. As it will create a json dictionairy.
So a .json file with as key the sequence id and as value either the seuence or taxonomy.

This is a lot more efficient than parsing through all the sequence ids again. 
For Bacteria and Archaea create a json per subset. 
