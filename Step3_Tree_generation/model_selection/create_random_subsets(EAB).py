import json
from Bio import SeqIO
import random
#if from fasta file directly:
CPA_id_list= []
#to change it, just edit the hmmaligned sequences and the name of the output file
HMMaligned = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/Bacteria_passed_all_filters_OKT3_GTDBreps_alignedPF00999.fasta'
outdir ='/home/mad149/00_nesi_projects/uc04105_nobackup/Likelihood_CPA'



for record in SeqIO.parse(HMMaligned, 'fasta'):
    CPA_id_list.append(record.id)
CPA_id_list = list(set(CPA_id_list))

def return_fasta_from_list(HMMaligned, random_list, count, outdir):
    written_list = []
    with open('{}/Bac_HMMscanned_PF00999_aligned_randomset_{}.fasta'.format(outdir, count), 'a') as OUT:
        for record in SeqIO.parse(HMMaligned, 'fasta'):
            if record.id not in random_list:
                pass
            else:
                if record.id not in written_list:
                    out_fasta = '>' + record.id + '\n' + str(record.seq) + '\n'
                    written_list.append(record.id)
                    OUT.write(out_fasta)
                    print('written {}'.format(count))
                else:
                    #I exclude entries as not to duplicate
                    print('nah')
                    pass
random_list_dict = {}

#pick and generate a random list

#when chosing 100 sequences out of 134, and order does not matter and no repetitions are made
# the total a
for count in range(1,101):
    #here we generate a random list
    random_list = random.sample(CPA_id_list, k=100)
    #we check if we did not accidentaly sample the same entry twice
    print(len(list(set(random_list))))
    if len(list(set(random_list))) != 100:
        #if we do have duplicates make a new random list and repeat the same procedure
        for i in random.sample(CPA_id_list, k=len(CPA_id_list)):
            if len(list(set(random_list))) < 100:
                if i not in random_list:
                    random_list.append(i)
                        
                else:
                    pass

    #if 100 entries are unqiue start writing and checking if we have made a similar list
    else:
        print('test')
        return_fasta_from_list(HMMaligned, random_list, count, outdir)
        pass
    
