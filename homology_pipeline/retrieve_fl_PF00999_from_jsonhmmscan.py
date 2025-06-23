import json
from Bio import SeqIO


scanned ='/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/Eukarya_PF00999_aligned_hmmscanned_onlyPF00999.json'
fl_fasta_in='/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/All_merged_pf00999_fulllength.faa'
fl_fasta_out='/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/All_merged_pf00999_fulllength_passed_scan.faa'

f = open(scanned,)
data=json.load(f)


match_list = []
written_list = []

with open(fl_fasta_out, 'w') as O:
    for record in SeqIO.parse(fl_fasta_in, 'fasta'):
        if record.id not in written_list:
            
            if record.id in list(data.keys()):
                fasta = '>' + record.id + '\n' + str(record.seq) + '\n' 
                O.write(fasta)
                written_list.append(record.id)
            else:
                pass
        else:
            pass
            
