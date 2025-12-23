import sys
import os
Arc_sequences_passed_thresholds = {}

from Bio import SeqIO

for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/Archaea_passed_all_filters_7OKT_GTDBreps_alignedPF00999.fasta','fasta'):
    Arc_sequences_passed_thresholds[str(record.id)] = record.seq
Hmmsearch = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Archaea'

sequences_hit_by_pf009999 = []

for hits in os.listdir(Hmmsearch):
    if 'PF00999' in hits:
        

        with open('{}/{}'.format(Hmmsearch, hits), 'r') as PF00999:
            for line in PF00999:
                if line[0] == '#':
                    pass
                else:
                    sequences_hit_by_pf009999.append(line.split(' ')[0])
        
                    
        sequences_hit_by_pf009999.extend(list(set(sequences_hit_by_pf009999)))
print(len(sequences_hit_by_pf009999))
    
    
list_pf00999_passed = []
for sequence in sequences_hit_by_pf009999:
    if sequence in  Arc_sequences_passed_thresholds:
        list_pf00999_passed.append(sequence)
print(len(list_pf00999_passed))

Hmmsearch = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Archaea'

sequences_hit_by_Bac = []

for hits in os.listdir(Hmmsearch):
    if 'Bac' in hits:
        

        with open('{}/{}'.format(Hmmsearch, hits), 'r') as Bac:
            for line in Bac:
                if line[0] == '#':
                    pass
                else:
                    sequences_hit_by_Bac.append(line.split(' ')[0])
        
                    
        sequences_hit_by_Bac.extend(list(set(sequences_hit_by_Bac)))
print(len(sequences_hit_by_Bac))
    
    
list_Bac_passed = []
for sequence in sequences_hit_by_Bac:
    if sequence in  Arc_sequences_passed_thresholds:
        list_Bac_passed.append(sequence)
print(len(list_Bac_passed))

Hmmsearch = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Archaea'

sequences_hit_by_Arc = []

for hits in os.listdir(Hmmsearch):
    if 'ArcHMM' in hits:
        

        with open('{}/{}'.format(Hmmsearch, hits), 'r') as Arc:
            for line in Arc:
                if line[0] == '#':
                    pass
                else:
                    sequences_hit_by_Arc.append(line.split(' ')[0])
        
                    
        sequences_hit_by_Arc.extend(list(set(sequences_hit_by_Arc)))
print(len(sequences_hit_by_Arc))
    
    
list_Arc_passed = []
for sequence in sequences_hit_by_Arc:
    if sequence in  Arc_sequences_passed_thresholds:
        list_Arc_passed.append(sequence)
print(len(list_Arc_passed))

Hmmsearch = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Archaea'

sequences_hit_by_Euk = []

for hits in os.listdir(Hmmsearch):
    if 'Euk' in hits:
        with open('{}/{}'.format(Hmmsearch, hits), 'r') as Euk:
            for line in Euk:
                if line[0] == '#':
                    pass
                else:
                    sequences_hit_by_Euk.append(line.split(' ')[0])
        
                    
        sequences_hit_by_Euk.extend(list(set(sequences_hit_by_Euk)))
print(len(sequences_hit_by_Euk))
    
    
list_Euk_passed = []
for sequence in sequences_hit_by_Euk:
    if sequence in Arc_sequences_passed_thresholds:
        list_Euk_passed.append(sequence)
print(len(list_Euk_passed))

#misses taxonomy bit
