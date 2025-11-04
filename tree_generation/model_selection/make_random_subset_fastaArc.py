import json
from Bio import SeqIO
import random
#if from fasta file directly:
CPA_id_list= []
HMMaligned = '/nesi/nobackup/uc04105/fasta_files/Seed_NhaA/SG_seed_update_jan142024_135seq_PF00999aligned.fasta.fasta'
outdir ='/nesi/nobackup/uc04105/fasta_files/Seed_NhaA'
for record in SeqIO.parse(HMMaligned, 'fasta'):
   CPA_id_list.append(record.id.split('|')[1])
CPA_id_list = list(set(CPA_id_list))

def return_fasta_from_list(HMMaligned, random_list2, count, outdir):
    written_list = []
    with open('{}/Seed_HMMscanned_PF00999_aligned_randomset_{}.fasta'.format(outdir, count), 'a') as OUT:
        for record in SeqIO.parse(HMMaligned, 'fasta'):
            if record.id.split('|')[1] not in random_list2:
                pass
            else:
                if record.id.split('|')[1] not in written_list:
                    out_fasta = '>' + record.id.split('|')[1] + '\n' + str(record.seq) + '\n'
                    written_list.append(record.id.split('|')[1])
                    OUT.write(out_fasta)
                else:
                    pass

for count in range(1,101):
    random_list = random.choices(CPA_id_list, k=100)
    if len(list(set(random_list))) != 100:
        for i in random.sample(CPA_id_list, k=len(CPA_id_list)):
            if len(list(set(random_list))) < 100:
                random_list.append(i)
            elif len(list(set(random_list))) == 100:
                break
    random_list2 = list(set(random_list))
    print(len(random_list2))
    random_list = []
    return_fasta_from_list(HMMaligned, random_list2, count, outdir)
    
    return_fasta_from_list(HMMaligned, random_list, count)
    
