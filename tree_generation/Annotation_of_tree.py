from Bio import SeqIO
import json
import os
import random


#MMseq
from collections import Counter
MMseq_dir = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq'
ARCMMseq ='/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Archaea_allhmmscanned_final_clustered_at0.7.fasta_rep_seq.fasta'
BACMMseq ='/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Bacteria_allhmmscanned_final_clustered_at0.7.fasta_rep_seq.fasta'
EukMMseq ='/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Euk_allhmmscanned_final_clustered_at0.7.fasta_rep_seq.fasta'

Euk_cluster_list =[]
with open(EukMMseq, 'r') as E:
    for line in E:
        Euk_cluster_list.append(line.split('\t')[0])

Euk_cluster = Counter(Euk_cluster_list)

Arc_cluster_list =[]
with open(ARCMMseq, 'r') as E:
    for line in E:
        Arc_cluster_list.append(line.split('\t')[0])
Arc_cluster = Counter(Arc_cluster_list)



Bac_cluster_list =[]
with open(BACMMseq, 'r') as E:
    for line in E:
        Bac_cluster_list.append(line.split('\t')[0])
Bac_cluster = Counter(Bac_cluster_list)

#HMMalign

#NhaA
#TRK_N
#TRK_C
#CHX17
#CHX17_C
#cNMP
#Nha1_C
#USP

#Bacteria
def Parse_hmmalign_to_list(infile):
    dom_dic =[]
    for record in SeqIO.parse('{}'.format(infile), 'fasta'):
        dom_dic.append(record.id)
    return(dom_dic)


#NhaA
Bac_NhaA = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/06965_vsAllBacteria_06965_aligned.fa.fasta')
Arc_NhaA = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/06965_vsAllArchaea_06965_aligned.fa.fasta')
Euk_NhaA = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/06965_vsAllEukarya_06965_aligned.fa.fasta')

#TRK_N
Bac_TRKN = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/02254_vsAllBacteria_02254_aligned.fa.fasta')
Arc_TRKN = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/02254_vsAllArchaea_02254_aligned.fa.fasta')
Euk_TRKN = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/02254_vsAllEukarya_02254_aligned.fa.fasta')

#TRK_C
Bac_TRKC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/02080_vsAllBacteria_02080_aligned.fa.fasta')
Arc_TRKC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/02080_vsAllArchaea_02080_aligned.fa.fasta')
Euk_TRKC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/02080_vsAllEukarya_02080_aligned.fa.fasta')

#CHX17
Bac_CH17X = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/23256_vsAllBacteria_23256_aligned.fa.fasta')
Arc_CH17X = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/23256_vsAllArchaea_23256_aligned.fa.fasta')
Euk_CH17X = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/23256_vsAllEukarya_23256_aligned.fa.fasta')

#CHX17_C
Bac_CH17XC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/23259_vsAllBacteria_23259_aligned.fa.fasta')
Arc_CH17XC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/23259_vsAllArchaea_23259_aligned.fa.fasta')
Euk_CH17XC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/23259_vsAllEukarya_23259_aligned.fa.fasta')


#Nha1_C
Bac_Nha1_C = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/08619_vsAllBacteria_08619_aligned.fa.fasta')
Arc_Nha1_C = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/08619_vsAllArchaea_08619_aligned.fa.fasta')
Euk_Nha1_C = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/08619_vsAllEukarya_08619_aligned.fa.fasta')


#cNMP
Bac_cNMP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/00027_vsAllBacteria_00027_aligned.fa.fasta')
Arc_cNMP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/00027_vsAllArchaea_00027_aligned.fa.fasta')
Euk_cNMP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/00027_vsAllEukarya_00027_aligned.fa.fasta')


#USP
Bac_USP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/00582_vsAllBacteria_00582_aligned.fa.fasta')
Arc_USP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/00582_vsAllArchaea_00582_aligned.fa.fasta')
Euk_USP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/00582_vsAllEukarya_00582_aligned.fa.fasta')

#tax_list = []
#print(len(Bac_USP))
#for i in Bac_USP:
#    tax_list.append(i.split('tax:')[1].split(';c_')[0])
#print(Counter(tax_list))

import os

Prot_dict = {}
count = 0
hmmaligned = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/cross_domain_final'
for aligned_set in os.listdir(hmmaligned):
    if 'faa.fasta'  in aligned_set:
        for record in SeqIO.parse('{}/{}'.format(hmmaligned, aligned_set), 'fasta'):
            seq = str(record.seq)
            tax = record.id.split('tax:')[1]
            prot_id =  record.id.split('_tax')[0]
            Prot_dict[prot_id] = tax, seq
hmmaligned = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain/'
for aligned_set in os.listdir(hmmaligned):
    if 'faa.fasta'  in aligned_set:
        for record in SeqIO.parse('{}/{}'.format(hmmaligned, aligned_set), 'fasta'):
            seq = str(record.seq)
            tax = record.id.split('tax:')[1]
            prot_id =  record.id.split('_tax')[0]
            Prot_dict[prot_id] = tax, seq
hmmaligned = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain'
for aligned_set in os.listdir(hmmaligned):
    if 'fa.fasta'  in aligned_set:
        for record in SeqIO.parse('{}/{}'.format(hmmaligned, aligned_set), 'fasta'):
            seq = str(record.seq)
            tax = record.id.split('tax:')[1]
            prot_id =  record.id.split('_tax')[0]
            Prot_dict[prot_id] = tax, seq

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_full_tree_aligned_annotation_25JULI.tsv', 'a') as C:

    C.write('Protein_id' + '\t' + 'Sequence' + '\t' + 'ngaps' + '\t' + 'Motif' + '\t' + 'Domain' + '\t' + 'GTDB_taxonomy' + '\t' + 'Phyla' + '\t' + 'Seed' + '\t' + 'First_site' + '\t' + 'Second_site' + '\t' + 'last_site' + '\t' + 'Cluster_count' +  '\t' + 'USP' + '\t' + 'NhaA' + '\t' + 'TRK_N' + '\t' + 'TRK_C' + '\t' + 'CHX17' + '\t' + 'CHX17_C' + '\t' + 'cNMP' + '\t' + 'Nha1_C'+ '\n')

    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/SG_seed_update_jan142024_135seq_PF00999aligned.fasta.fasta', 'fasta'):
        seq = str(record.seq)
        tax = record.id.split('|')[1]
        prot_id =  record.id.split('_')[-1] + '_' + record.id.split('_')[-2] + '_' + str(count)
        count = count + 1
        Prot_dict[prot_id] = tax, seq

    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/tree_seqs/CPA_Phylogney_hmmscanned_hmmaligned_PF00999_25Juli.faa', 'fasta'):
        NhaA = 'NA'
        TRK_N = 'NA'
        TRK_C = 'NA'
        CHX17 = 'NA'
        CHX17_C = 'NA'
        cNMP = 'NA'
        Nha1_C = 'NA'
        USP = 'NA'
        if (record.id.split('_')[0]) == 'Arc226':
            Domain = 'Archaea'

            Seed = ''
            GTDB_taxonomy = Prot_dict[record.id][0]
            Phyla = GTDB_taxonomy.split(';')[1]
            Class = GTDB_taxonomy.split(';')[2]
            Order = GTDB_taxonomy.split(';')[3]
            Family = GTDB_taxonomy.split(';')[4]
            Genus = GTDB_taxonomy.split(';')[5]
            Species = GTDB_taxonomy.split(';s__')[1]
            Cluster_count = Arc_cluster[(record.id+'_tax:{}'.format(GTDB_taxonomy))]

            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_NhaA:
                NhaA = 'NhaA'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_TRKN:
                TRK_N = 'TRK_N'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_TRKC:
                TRK_C = 'TRK_C'
                combination_TRK = TRK_N + TRK_C
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_CH17X:
                CHX17 = 'CHX17'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_CH17XC:
                CHX17_C = 'CHX17_C'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_cNMP:
                cNMP = 'cNMP'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_Nha1_C:
                Nha1_C = 'Nha1_C'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Arc_USP:
                USP = 'USP'
                
        elif (record.id.split('_')[0]) == 'Bac226':
            Domain = 'Bacteria'
            Seed = ''
            GTDB_taxonomy = Prot_dict[record.id][0]
            Phyla = GTDB_taxonomy.split(';')[1]
            Class = GTDB_taxonomy.split(';')[2]
            Order = GTDB_taxonomy.split(';')[3]
            Family = GTDB_taxonomy.split(';')[4]
            Genus = GTDB_taxonomy.split(';')[5]
            Species = GTDB_taxonomy.split(';s__')[1]
            Cluster_count = Bac_cluster[(record.id+'_tax:{}'.format(GTDB_taxonomy))]

            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_NhaA:
                NhaA = 'NhaA'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_TRKN:
                TRK_N = 'TRK_N'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_TRKC:
                TRK_C = 'TRK_C'
                combination_TRK = TRK_N + TRK_C
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_CH17X:
                CHX17 = 'CHX17'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_CH17XC:
                CHX17_C = 'CHX17_C'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_cNMP:
                cNMP = 'cNMP'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_Nha1_C:
                Nha1_C = 'Nha1_C'
            if record.id+'_tax:{}'.format(GTDB_taxonomy) in Bac_USP:
                USP = 'USP'
        elif record.id.split('_')[0] == 'EukM6':
            Domain = 'Eukarya'
            Seed = ''
            GTDB_taxonomy = Prot_dict[record.id][0]
            Phyla = ''
            Class = ''
            Order =''
            Family =''
            Genus = ''
            Species = GTDB_taxonomy.split(';')[-1]

            Cluster_count = Euk_cluster[(record.id+'_tax:{}'.format(Species))]

            if record.id+'_tax:{}'.format(Species) in Euk_NhaA:
                NhaA = 'NhaA'
            if record.id+'_tax:{}'.format(Species) in Euk_TRKN:
                TRK_N = 'TRK_N'
            if record.id+'_tax:{}'.format(Species) in Euk_TRKC:
                TRK_C = 'TRK_C'
                combination_TRK = TRK_N + TRK_C
            if record.id+'_tax:{}'.format(Species) in Euk_CH17X:
                CHX17 = 'CHX17'
            if record.id+'_tax:{}'.format(Species) in Euk_CH17XC:
                CHX17_C = 'CHX17_C'
            if record.id+'_tax:{}'.format(Species) in Euk_cNMP:
                cNMP = 'cNMP'
            if record.id+'_tax:{}'.format(Species) in Euk_Nha1_C:
                Nha1_C = 'Nha1_C'
            if record.id+'_tax:{}'.format(Species) in Euk_USP:
                USP = 'USP'
            
            
        else:
            Domain = 'Seed'
            Seed = record.id
            GTDB_taxonomy = Prot_dict[record.id][0]
            Phyla = ''
            Class = ''
            Order =''
            Family =''
            Genus = ''
            Species = GTDB_taxonomy.split(';')[-1]
            Cluster_count = 'NA'
            NhaA = 'NA'
            TRK_N = 'NA'
            TRK_C = 'NA'
            CHX17 = 'NA'
            CHX17_C = 'NA'
            cNMP = 'NA'
            Nha1_C = 'NA'
            USP = 'NA'
            
        ngaps = str(record.seq.count('-'))
        First_site = (record.seq[145])
        Second_site = record.seq[146]
        last_site = record.seq[300]
        Motif = record.seq[145] + record.seq[146] + record.seq[300]
        C.write(record.id + '\t' + str(record.seq) + '\t' + ngaps + '\t' + Motif + '\t' + Domain + '\t' + GTDB_taxonomy + '\t' + Phyla + '\t' + Seed + '\t' + First_site + '\t' + Second_site + '\t' + last_site + '\t' + str(Cluster_count) +  '\t' + USP + '\t' + NhaA + '\t' + TRK_N + '\t' + TRK_C + '\t' + CHX17 + '\t' + CHX17_C + '\t' + cNMP + '\t' + Nha1_C + '\n')

