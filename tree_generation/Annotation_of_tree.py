from Bio import SeqIO
import json
import os
import random


#MMseq

#count how many sequences a cluster rep actually represents.
from collections import Counter
ARCMMseq ='/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Archaea_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_cluster.tsv'
BACMMseq ='/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Bacteria_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_cluster.tsv'
EukMMseq ='/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Eukarya_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_cluster.tsv'

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


#Bacteria
def Parse_hmmalign_to_list(infile):
    dom_dic =[]
    for record in SeqIO.parse('{}'.format(infile), 'fasta'):
        if '|' not in record.id:
            dom_dic.append(record.id)
        else:
            dom_dic.append(record.id.split('|')[1])
    return(dom_dic)

#NhaA
Bac_NhaA = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF06965_NhaA.fasta')
Arc_NhaA = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF06965_NhaA.fasta')
Euk_NhaA = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF06965_NhaA.fasta')

#TRK_N
Bac_TRKN = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF02254_TrkN.fasta')
Arc_TRKN = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF02254_TrkN.fasta')
Euk_TRKN = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF02254_TrkN.fasta')

#TRK_C
Bac_TRKC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF02080_TrkC.fasta')
Arc_TRKC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF02080_TrkC.fasta')
Euk_TRKC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF02080_TrkC.fasta')

#CHX17
Bac_CH17X = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF23256_CHX2nd.fasta')
Arc_CH17X = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF23256_CHX2nd.fasta')
Euk_CH17X = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF23256_CHX2nd.fasta')

#CHX17_C
Bac_CH17XC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF23259_CHXC.fasta')
Arc_CH17XC = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF23259_CHXC.fasta')
Euk_CH17XC =Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF23259_CHXC.fasta')
#Nha1_C
Bac_Nha1_C = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF08619_NhaC1term.fasta')
Arc_Nha1_C = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF08619_NhaC1term.fasta')
Euk_Nha1_C = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF08619_NhaC1term.fasta')
#cNMP
Bac_cNMP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Bacteria_fl_vsPF00027_cNMP.fasta')
Arc_cNMP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Archaea_fl_vsPF00027_cNMP.fasta')
Euk_cNMP = Parse_hmmalign_to_list('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/HMMalign/Eukarya_fl_vsPF00027_cNMP.fasta')
#Use last big file with taxonomic annotation

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


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_full_tree_aligned_annotation_29thOct.tsv', 'w') as C:

    C.write('Protein_id' + '\t' + 'Sequence' + '\t' + 'ngaps' + '\t' + 'Motif' + '\t' + 'Domain' + '\t' + 'GTDB_taxonomy' + '\t' + 'Phyla' + '\t' + 'Seed' + '\t' + 'First_site' + '\t' + 'Second_site' + '\t' + 'last_site' + '\t' + 'Cluster_count' '\t' + 'NhaA' + '\t' + 'TRK_N' + '\t' + 'TRK_C' + '\t' + 'CHX17' + '\t' + 'CHX17_C' + '\t' + 'cNMP' + '\t' + 'Nha1_C'+ '\n')

    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/SG_seed_update_jan142024_135seq_alignedPF009999_filtered.faa.fasta', 'fasta'):
        seq = str(record.seq)
        if '_d' in record.id:
            tax = record.id.split('_d')[1]
        elif '___' in record.id:
            tax = record.id.split('___')[1]
        else:
            tax = 'Eukarya'
        prot_id =  record.id.split('_')[-1] + '_' + record.id.split('_')[-2] 
        count = count + 1
        Prot_dict[prot_id] = tax, seq

    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/CPA_Phylogeny_redone_6Okt_PF00999_aligned_hmmscanned.fasta', 'fasta'):
        NhaA = 'NA'
        TRK_N = 'NA'
        TRK_C = 'NA'
        CHX17 = 'NA'
        CHX17_C = 'NA'
        cNMP = 'NA'
        Nha1_C = 'NA'
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
            Cluster_count = Arc_cluster[record.id]

            if record.id in Arc_NhaA:
                NhaA = 'NhaA'
            if record.id in Arc_TRKN:
                TRK_N = 'TRK_N'
            if record.id in Arc_TRKC:
                TRK_C = 'TRK_C'
                combination_TRK = TRK_N + TRK_C
            if record.id in Arc_CH17X:
                CHX17 = 'CHX17'
            if record.id in Arc_CH17XC:
                CHX17_C = 'CHX17_C'
            if record.id in Arc_cNMP:
                cNMP = 'cNMP'
            if record.id in Arc_Nha1_C:
                Nha1_C = 'Nha1_C'

                
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
            Cluster_count = Bac_cluster[record.id]

            if record.id in Bac_NhaA:
                NhaA = 'NhaA'
            if record.id in Bac_TRKN:
                TRK_N = 'TRK_N'
            if record.id in Bac_TRKC:
                TRK_C = 'TRK_C'
                combination_TRK = TRK_N + TRK_C
            if record.id in Bac_CH17X:
                CHX17 = 'CHX17'
            if record.id in Bac_CH17XC:
                CHX17_C = 'CHX17_C'
            if record.id in Bac_cNMP:
                cNMP = 'cNMP'
            if record.id in Bac_Nha1_C:
                Nha1_C = 'Nha1_C'

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

            Cluster_count = Euk_cluster[record.id]

            if record.id in Euk_NhaA:
                NhaA = 'NhaA'
            if record.id in Euk_TRKN:
                TRK_N = 'TRK_N'
            if record.id in Euk_TRKC:
                TRK_C = 'TRK_C'
                combination_TRK = TRK_N + TRK_C
            if record.id in Euk_CH17X:
                CHX17 = 'CHX17'
            if record.id in Euk_CH17XC:
                CHX17_C = 'CHX17_C'
            if record.id in Euk_cNMP:
                cNMP = 'cNMP'
            if record.id in Euk_Nha1_C:
                Nha1_C = 'Nha1_C'

            
            
        else:
            Domain = 'Seed'
            Seed = record.id
            GTDB_taxonomy = Prot_dict[record.id.split('_')[-1] + '_' + record.id.split('_')[-2]][0]
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

            
        ngaps = str(record.seq.count('-'))
        First_site = (record.seq[145])
        Second_site = record.seq[146]
        last_site = record.seq[300]
        Motif = record.seq[145] + record.seq[146] + record.seq[300]
        C.write(record.id + '\t' + str(record.seq) + '\t' + ngaps + '\t' + Motif + '\t' + Domain + '\t' + GTDB_taxonomy + '\t' + Phyla + '\t' + Seed + '\t' + First_site + '\t' + Second_site + '\t' + last_site + '\t' + str(Cluster_count) +  '\t' + NhaA + '\t' + TRK_N + '\t' + TRK_C + '\t' + CHX17 + '\t' + CHX17_C + '\t' + cNMP + '\t' + Nha1_C + '\n')

#check if annotation and alignment match

from Bio import SeqIO

alignment_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/CPA_Phylogeny_redone_6Okt_PF00999_aligned_hmmscanned.fasta', 'fasta'):
    alignment_list.append(record.id)

print(len(alignment_list))

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_full_tree_aligned_annotation_29thOct.tsv', 'r') as A:
    next(A, None)
    for line in A:
        if (line.split('\t')[0]) not in alignment_list:
            print(line.split('\t')[0])
