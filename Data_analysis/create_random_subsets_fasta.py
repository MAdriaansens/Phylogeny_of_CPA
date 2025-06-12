import random
import json
from Bio import SeqIO

HMMaligned = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/All_merged.faa'
outdir='/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/random_list'
direct ='/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan'
HMS00999='Eukarya_PF00999_aligned_hmmscanned_onlyPF00999.json'
CPA_id_list = []

with open('{}/{}'.format(direct, HMS00999)) as json_file:
    CPA_data=json.load(json_file)
    for key in CPA_data.keys():
        output = CPA_data[key]
        if output[0] == 'PF00999':
            CPA_id_list.append(key.split('_tax:')[0])


print(len(list(set(CPA_id_list)))))
print(len(CPA_id_list))

def return_fasta_from_list(HMMaligned, random_list, count):
    with open('{}/Euk_HMMscanned_PF00999_aligned_randomset_{}.fasta'.format(outdir, count), 'a') as OUT:
        for record in SeqIO.parse(HMMaligned, 'fasta'):
            if record.id.split('_tax:')[0] in random_list:
                out_fasta = '>' + record.id + '\n' + str(record.seq) + '\n'
                print(out_fasta)
                OUT.write(out_fasta)

for count in range(1,101):
    random_list = random.choices(CPA_id_list, k=100)
    return_fasta_from_list(HMMaligned, random_list, count)
    
