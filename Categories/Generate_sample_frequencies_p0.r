library(dplyr)
library(tidyr)

IT_Bacteria <- read.delim('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/IT_Bacteria_9SEPT.tsv', sep='\t', stringsAsFactors = TRUE)

Frequencies <- IT_Bacteria %>% count(Sample)
write.table(Frequencies, file='Bacteria_sample_frequencies.tsv', quote=FALSE, sep='\t', col.names = NA)



IT_Archaea <- read.delim('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/IT_Archaea_8SEPT.tsv', sep='\t', stringsAsFactors = TRUE)

Frequencies <- IT_Archaea %>% count(Sample)
write.table(Frequencies, file='Archaea_sample_frequencies.tsv', quote=FALSE, sep='\t', col.names = NA)

#I then .lowercase the samples and only take the unique values. 
#I then manually annotate 9000 sample sites. ~800 in Archaea and 7800 in Bacteria
