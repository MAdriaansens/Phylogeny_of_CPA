import json
from Bio import SeqIO
import random
direct ='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea'
HMS00999='Archaea_PF00999_aligned_hmmscanned_onlyPF00999.json'
outdir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/intermediate_PF00999/random_list'


CPA_id_list=[]
with open('{}/{}'.format(direct, HMS00999)) as json_file:
    CPA_data=json.load(json_file)
    for key in CPA_data.keys():
        output = CPA_data[key]
        if output[0] == 'PF00999':
            CPA_id_list.append(key.split('_tax:')[0])

HMMaligned = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/intermediate_PF00999/PF00999_Archaea_MERGED_alignedPF00999.faa'

def return_fasta_from_list(HMMaligned, random_list, count):
    written_list = []
    with open('{}/Arc_HMMscanned_PF00999_aligned_randomset_{}.fasta'.format(outdir, count), 'a') as OUT:
        for record in SeqIO.parse(HMMaligned, 'fasta'):

            if record.id.split('_tax:')[0].split('|')[1] not in random_list:
                pass
            else:
                if record.id.split('_tax:')[0].split('|')[1] not in written_list:
                    out_fasta = '>' + record.id.split('|')[1] + '\n' + str(record.seq) + '\n'

                    written_list.append(record.id.split('_tax:')[0].split('|')[1])
                    OUT.write(out_fasta)
                else:
                    pass
for count in range(1,101):
    random_list = []
    random_list = random.choices(CPA_id_list, k=100)
    print(len(random_list))
    return_fasta_from_list(HMMaligned, random_list, count)
    
