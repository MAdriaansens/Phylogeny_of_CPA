output_file='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/Archaea_all_part1A_remove_dupes.fasta'
from Bio import SeqIO
input_file='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/Archaea_all_part1A.fasta'
Entry_list = []
Entry_dict = {}

for record in SeqIO.parse(input_file, 'fasta'):
    if record.id in Entry_list:
        pass
    else:
        Entry_list.append(record.id)
        Entry_dict[record.id] = str(record.seq)

with open(output_file, 'w') as output:
    for entry in set(Entry_list):
        line = '>' + entry + '\n' + Entry_dict[entry] + '\n'
        output.write(line)
