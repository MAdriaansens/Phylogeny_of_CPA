from Bio import SeqIO
import json
import os
Euk_meta = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv'

file_path = '/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax.json'

protein_ids_all = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Eukarya_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_all_seqs.fasta', 'fasta'):
    protein_ids_all.append(record.id)
if os.path.exists(file_path):
    print('yah')
    
    with open('/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax.json', 'r') as F:
       data = json.load(F)
else:
    data = {}
    #make a dictionairy for representatives and what their GTDB_id is
    with open(Euk_meta, 'r') as Euks:
        for Euk in Euks:
            protein_id = Euk.split('\t')[0]
            if protein_id in protein_ids_all:
                Tax_inf = Euk.split('\t')[1] +'_Taxa:' + Euk.split('\t')[4]
                data[protein_id] = Tax_inf 
    with open('/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax_tax.json', 'w') as J:
        json.dump(data, J)





MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Eukarya_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_cluster.tsv'

rep_list = []


with open(MMseq, 'r') as clusters:
    for cluster in clusters:
        rep = cluster.split('\t')[0]
        rep_list.append(rep)
    clusters.close()

from collections import Counter

#the rep list, is just a list of all values in the representative columns, this means that ids are present multiple times

rep_list = (rep_list)
reps = Counter(rep_list)

#reps returns is a counter of how many sequences a rep represents
items = reps.items()


singleton_list = []
multiple_list = []

for i in items:
    if i[-1] == 1:
        singleton_list.append(i[0])
    else:
        multiple_list.append(i[0])

print(len(multiple_list))
print(len(singleton_list))
replist = []
Rep_dict = {}

with open(MMseq, 'r') as clustermap:
    count = 0
    second_list = []
    element_count = 0
    for entry in clustermap:
        count = reps[entry.split('\t')[0]]
        repm = entry.split('\t')[0]

        if repm in multiple_list:
            if repm in replist:
                second = (entry.split('\t')[1].split('\n')[0])
                second_list.append(second)
                if len(second_list) == count:
                    Rep_dict[repm] = second_list
            else:

                replist.append(repm)
                second_list = []
                second_list.append(repm)
        else:
            pass
print(len(set(replist)))

 
Tax = []
count = 0
same_species_count = 0
same_species_count1 = 0
same_species_count2 = 0
same_species_count3 = 0
same_species_count4 = 0
same_species_count5 = 0
same_species_count6 = 0
same_species_count7 = 0
same_species_count8 = 0
same_species_count9=0
same_species_count10= 0
same_species_count11=0
same_species_count12= 0
same_species_count13=0
same_species_count14=0
same_species_count15=0
same_species_count16=0
same_species_count17=0
same_species_count18=0
num_same_species_proteins = 0
Taxonomy = ''
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Eukarya_clusters_7OKT.tsv', 'w') as Out:
    Header = 'Representative_id' + '\t' + 'No_proteins_in_cluster' + '\t' + 'No_diff_species' +  '\t' + 'Last_tax_grouping' +  '\t' + 'Major_tax' + '\t' +  'Shared_Taxonomy'  + '\t' + 'Taxa' +'\t' + 'Species_list' + '\t' + 'clusteroid_list' + '\n'
    Out.write(Header)
    for euk in replist:
        Shared_tax = ''
        minor = []
        minor_1 = []
        minor_2 = []
        minor_3 = []
        minor_4 =[]
        minor_5=[]
        minor_6=[]
        minor_7 = []
        minor_8=[]
        minor_9= []
        minor_10=[]
        minor_11= []
        minor_12=[]
        minor_13=[]
        minor_14=[]
        minor_15=[]
        minor_16=[]
        minor_17=[]
        Taxa_in_cluster_list = []
        clusteroid_list = Rep_dict[euk]
        Species_list = []
        Taxa = ''
        MJ = []

        for prot in clusteroid_list:
            Taxa_in_cluster_list.append(data[prot])
        for x in Taxa_in_cluster_list:
            MJ.append(x.split(';')[4])
           
        if len(set(MJ)) > 1:
            print(set(MJ))
            print(len(Taxa_in_cluster_list))
        MJ = str(list(set(MJ)))
        if len(set(Taxa_in_cluster_list)) == 1:
            Shared_tax = 'same_species'
            Representative_id = euk
            Taxonomy = Taxa_in_cluster_list[0]
            No_proteins_in_cluster = (len(clusteroid_list))
            List_proteins = clusteroid_list
            No_diff_species = (len(set(Taxa_in_cluster_list)))
            Taxa = 'equal as Taxonomy'
            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
            same_species_count +=1
        else:
            for i in Taxa_in_cluster_list:
                Species_list.append(i.split('_Taxa')[0])
                minor.append(i.split(';')[4])
                if i.split(';')[5] == 'Fungiincertaesedis':
                    minor_1.append('Fungi_incertae_sedis')
                    Tax.append(i.split(';')[4])
                        
                else:
                    minor_1.append(i.split(';')[5])
                    Tax.append(i.split(';')[4])
                minor_2.append(i.split(';')[6])
                minor_3.append(i.split(';')[7])
                minor_4.append(i.split(';')[8])
                if len(i.split(';')) > 9:
                
                    minor_5.append(i.split(';')[9])
                
                if len(i.split(';')) > 10:
                    minor_6.append(i.split(';')[10])
                if len(i.split(';')) > 11:
                    minor_7.append(i.split(';')[11])
                if len(i.split(';')) > 12:
                    minor_8.append(i.split(';')[12])
                if len(i.split(';')) > 13:
                    minor_9.append(i.split(';')[13])
                if len(i.split(';')) > 14:
                    minor_10.append(i.split(';')[14])
                if len(i.split(';')) > 15:
                    minor_11.append(i.split(';')[15])
                if len(i.split(';')) > 16:
                    minor_12.append(i.split(';')[16])
                if len(i.split(';')) > 17:
                    minor_13.append(i.split(';')[17])
                if len(i.split(';')) > 18:
                    minor_14.append(i.split(';')[18])
                if len(i.split(';')) > 19:
                    minor_15.append(i.split(';')[19])
                if len(i.split(';')) > 20:
                    minor_16.append(i.split(';')[20])
                if len(i.split(';')) > 21:
                    minor_17.append(i.split(';')[21])
                    pass
            if len(set(minor)) != 1:           
                Shared_tax = format(i.split(';')[3])
                Representative_id = euk
                No_proteins_in_cluster = (len(clusteroid_list))
                List_proteins = clusteroid_list
                No_diff_species = (len(set(Taxa_in_cluster_list)))
                Taxonomy = 'Eukarya'
                Taxa =  str(set(minor))
                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                same_species_count1 +=1
                print(Line)
            else:
                if (len(set(minor_1))) != 1:
                    Shared_tax = (i.split(';')[4])
                    Representative_id = euk
                    No_proteins_in_cluster = (len(clusteroid_list))
                    List_proteins = clusteroid_list
                    No_diff_species = (len(set(Taxa_in_cluster_list)))
                    Taxonomy = 'Eukarya_' + minor[0]
                    Taxa =  str(set(minor_1))
                    Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                    same_species_count2 +=1
    
                else:
                    if (len(set(minor_2))) != 1:
                        Shared_tax = (i.split(';')[5])
                        Representative_id = euk
                        No_proteins_in_cluster = (len(clusteroid_list))
                        List_proteins = clusteroid_list
                        No_diff_species = (len(set(Taxa_in_cluster_list)))
                        Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0]
                        Taxa =  str(set(minor_2))
                        Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                        same_species_count3 +=1
                    else:
                        if len(set(minor_3)) !=1:
                            Shared_tax = (i.split(';')[6])
                            Representative_id = euk
                            No_proteins_in_cluster = (len(clusteroid_list))
                            List_proteins = clusteroid_list
                            No_diff_species = (len(set(Taxa_in_cluster_list)))
                            Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0]
                            Taxa =  str(set(minor_3))
                            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                            same_species_count4 +=1

                        else:
                            if len(set(minor_4)) !=1:
                                Shared_tax = (i.split(';')[7])
                                Representative_id = euk
                                No_proteins_in_cluster = (len(clusteroid_list))
                                List_proteins = clusteroid_list
                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0]
                                Taxa =  str(set(minor_4))
                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                same_species_count5 +=1
                            else:
                                if len(set(minor_5)) !=1:
                                    Shared_tax = (i.split(';')[8])
                                    Representative_id = euk
                                    No_proteins_in_cluster = (len(clusteroid_list))
                                    List_proteins = clusteroid_list
                                    Taxa =  str(set(minor_5))
                                    No_diff_species = (len(set(Taxa_in_cluster_list)))
                                    Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0]
                                    
                                    Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                    same_species_count6 +=1
                                else:
                                    if len(set(minor_6)) !=1:
                                        Shared_tax = (i.split(';')[9])
                                        Representative_id = euk
                                        No_proteins_in_cluster = (len(clusteroid_list))
                                        List_proteins = clusteroid_list
                                        Taxa =  str(set(minor_6))
                                        No_diff_species = (len(set(Taxa_in_cluster_list)))
                                        Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + minor_5[0]
                                        
                                        Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                        same_species_count7 +=1
                                    else:
                                         if len(set(minor_7)) !=1:
                                            Shared_tax = (i.split(';')[10])
                                            Representative_id = euk
                                            No_proteins_in_cluster = (len(clusteroid_list))
                                            List_proteins = clusteroid_list
                                            Taxa =  str(set(minor_7))
                                            No_diff_species = (len(set(Taxa_in_cluster_list)))
                                            Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0]
                                            
                                            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                            same_species_count8 +=1/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax_tax.json
                                         else:
                                             if len(set(minor_8)) !=1:
                                                Shared_tax = (i.split(';')[11])
                                                Representative_id = euk
                                                No_proteins_in_cluster = (len(clusteroid_list))
                                                List_proteins = clusteroid_list
                                                Taxa =  str(set(minor_8))
                                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0]
                                                
                                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                same_species_count9 +=1
                                             else:
                                                 if len(set(minor_9)) !=1:
                                                    Shared_tax = (i.split(';')[12])
                                                    Representative_id = euk
                                                    No_proteins_in_cluster = (len(clusteroid_list))
                                                    List_proteins = clusteroid_list
                                                    Taxa =  str(set(minor_9))
                                                    No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                    Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]
                                                    
                                                    Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                    same_species_count10 +=1
                                                    
                                                 else:
                                                    if len(set(minor_10)) !=1:
                                                        Shared_tax = (i.split(';')[13])
                                                        Representative_id = euk
                                                        No_proteins_in_cluster = (len(clusteroid_list))
                                                        List_proteins = clusteroid_list
                                                        Taxa =  str(set(minor_10))
                                                        No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                        Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0] + minor_9[0]
                                                        
                                                        Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                        same_species_count11 +=1
                                                    else:
                                                        if len(set(minor_11)) !=1:
                                                            Shared_tax = (i.split(';')[14])
                                                            Representative_id = euk
                                                            No_proteins_in_cluster = (len(clusteroid_list))
                                                            List_proteins = clusteroid_list
                                                            Taxa =  str(set(minor_11))
                                                            No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                            Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0] + '_' + minor_9[0]  + '_'+ minor_10[0]
                                                            
                                                            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                            same_species_count12 +=1
                                                        else:
                                                            if len(set(minor_12)) !=1:
                                                                Shared_tax = (i.split(';')[15])
                                                                Representative_id = euk
                                                                No_proteins_in_cluster = (len(clusteroid_list))
                                                                List_proteins = clusteroid_list
                                                                Taxa =  str(set(minor_12))
                                                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                                Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]+ '_' + minor_9[0] + '_'+ minor_10[0] + '_'+ minor_11[0]
                                                                
                                                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                                same_species_count13 +=1
                                                            else:
                                                                if len(set(minor_13)) !=1:
                                                                    Shared_tax = (i.split(';')[16])
                                                                    Representative_id = euk
                                                                    No_proteins_in_cluster = (len(clusteroid_list))
                                                                    List_proteins = clusteroid_list
                                                                    Taxa =  str(set(minor_13))
                                                                    No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                                    Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]+ '_' + minor_9[0] + '_'+ minor_10[0] + '_'+ minor_11[0] +'_' + minor_12[0]
                                                                    
                                                                    Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                                    same_species_count14 +=1
                                                                else:
                                                                    if len(set(minor_14)) !=1:
                                                                        Shared_tax = (i.split(';')[17])
                                                                        Representative_id = euk
                                                                        No_proteins_in_cluster = (len(clusteroid_list))
                                                                        List_proteins = clusteroid_list
                                                                        Taxa =  str(set(minor_14))
                                                                        No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                                        Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]+ '_' + minor_9[0] + '_'+ minor_10[0] + '_'+ minor_11[0] +'_' + minor_12[0] + '_' + minor_13[0]
                                                                        
                                                                        Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                                        same_species_count15 +=1
                                                                    else:
                                                                        if len(set(minor_15)) !=1:
                                                                            Shared_tax = (i.split(';')[18])
                                                                            Representative_id = euk
                                                                            No_proteins_in_cluster = (len(clusteroid_list))
                                                                            List_proteins = clusteroid_list
                                                                            Taxa =  str(set(minor_15))
                                                                            No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                                            Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]+ '_' + minor_9[0] + '_'+ minor_10[0] + '_'+ minor_11[0] +'_' + minor_12[0] + '_' + minor_13[0] + minor_14[0]
                                                                            
                                                                            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                                            same_species_count16 +=1
                                                                        else:
                                                                                
                                                                            if len(set(minor_16)) !=1:
                                                                                Shared_tax = (i.split(';')[19])
                                                                                Representative_id = euk
                                                                                No_proteins_in_cluster = (len(clusteroid_list))
                                                                                List_proteins = clusteroid_list
                                                                                Taxa =  str(set(minor_16))
                                                                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                                                Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]+ '_' + minor_9[0] + '_'+ minor_10[0] + '_'+ minor_11[0] +'_' + minor_12[0] + '_' + minor_13[0] + '_' + minor_14[0] + '_' + minor_15[0]
                                                                                
                                                                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                                                same_species_count17 +=1
                                                                            else:
                                                                                Shared_tax = (i.split(';')[20])
                                                                                Representative_id = euk
                                                                                No_proteins_in_cluster = (len(clusteroid_list))
                                                                                List_proteins = clusteroid_list
                                                                                Taxa =  str(set(minor_17))
                                                                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                                                Taxonomy =  'Eukarya_' + minor[0] + '_' + minor_1[0] +'_' + minor_2[0] + '_' + minor_3[0] + '_' + minor_4[0] + '_'  + minor_5[0] + '_' + minor_6[0] + '_' + minor_7[0] + '_' + minor_8[0]+ '_' + minor_9[0] + '_'+ minor_10[0] + '_'+ minor_11[0] +'_' + minor_12[0] + '_' + minor_13[0] + '_' + minor_14[0] + '_' + minor_15[0] + minor_16[0]
                                                                                
                                                                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(No_diff_species) +  '\t' + Shared_tax + '\t' + MJ + '\t' +  Taxonomy + '\t' + Taxa +'\t' + str(set(Species_list)) + '\t' + str(clusteroid_list) + '\n'
                                                                                same_species_count18 +=1
        Out.write(Line)
        count +=1
print(count)
