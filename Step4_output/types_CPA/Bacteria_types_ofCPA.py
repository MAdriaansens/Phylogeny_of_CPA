
MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_april28/tree_input/Bacteria_April_hmmaligned_e03_mmseq_treeinput_clusterd_at_0.7.faa_cluster.tsv'
rep_list = []


with open(MMseq, 'r') as clusters:
    for cluster in clusters:
        
        rep = cluster.split('\t')[0].split(':')[0]
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
print(len(items))

for i in items:
    if i[-1] == 1:
        singleton_list.append(i[0].split(':')[0])
    else:
        multiple_list.append(i[0].split(':')[0])

print(len(multiple_list))
print(len(singleton_list))
replist = []
Reps_dic = {}
with open(MMseq, 'r') as clustermap:
    count = 0
    total = 0
    second_list = []
    element_count = 0
    total = 0
    for entry in clustermap:
        total = total + 1
        count = reps[entry.split('\t')[0].split(':')[0]]
        repm = entry.split('\t')[0].split(':')[0]
        if repm in multiple_list:
            if repm in replist:
                second = (entry.split('\t')[1].split('\n')[0])
                second_list.append(second)
                if len(second_list) == count:
                    Reps_dic[repm] = second_list
            else:

                replist.append(repm)
                second_list = []
                second_list.append(repm)
        else:
            pass
print(total)
#load matching types
clade_dic = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_april28/Red_informed_clade_20_lg_cat_gamma_RED_interval_28April.tsv', 'r') as LG_Cat_Gamma_clade:
    next(LG_Cat_Gamma_clade, None)
    for line in LG_Cat_Gamma_clade:
        protein_id = line.split('\t')[0]
        domain = line.split('\t')[2]
        clade_dic[protein_id] = domain

seq_dic = {}
from Bio import SeqIO
tax_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bacteria_CPA_fl.fasta','fasta'):
    seq_dic[record.id.split(':')[0]] = [record.id, str(record.seq)]
    tax_list.append(record.id.split('tax:')[1])
print(len(tax_list))
tax_list = list(set(tax_list))
print(len(tax_list))











#list
CPA1_list = []
CHX_list = []
Kef_list = []
NhaA_list = []
Uncharacterized_list = []
NhaS5_list = []
Undescribed_CPA1_list = []
SOD2_list = []
NhaP_CPA1_list = []
UndProkarya_CPA1_IDK_list = []
GerN_list = []
CPA1_SL_list = []

for key in clade_dic.keys():

    if 'Bac226' in key:
        if key[0].isnumeric():
            pass
        else:

            if clade_dic[key] == '8':
                if key in Reps_dic:
                    for entry.split(':')[0] in Reps_dic[key]:
                        
                        
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_CPA1_sequences.faa", "a") as CPA1out:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            CPA1_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:CPA1'.format(id_tax) + '\n' + str(sequence) + '\n'
                            CPA1out.write(line)
                        CPA1out.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_CPA1_sequences.faa", "a") as CPA1out:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        CPA1_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:CPA1'.format(id_tax) + '\n' + str(sequence) + '\n'
                        CPA1out.write(line)
                    CPA1out.close()
            elif clade_dic[key] == '62243':
                    if key in Reps_dic:
                        for entry.split(':')[0] in Reps_dic[key]:
                            
                            with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_CHX_sequences.faa", "a") as CHXout:
                                sequence  = seq_dic[entry.split(':')[0]][1]
                                id_tax =  seq_dic[entry.split(':')[0]][0]
                                CHX_list.append(id_tax.split('tax:')[1])

                                line = '>{}_Protein:CHX'.format(id_tax) + '\n' + str(sequence) + '\n'
                                CHXout.write(line)
                            CHXout.close()
                    else:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_CHX_sequences.faa", "a") as CHXout:
                            sequence  = seq_dic[key][1]
                            id_tax =  seq_dic[key][0]
                            CHX_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:CHX'.format(id_tax) + '\n' + str(sequence) + '\n'
                            CHXout.write(line)
                        CHXout.close()
            elif clade_dic[key] == '29930':
                    if key in Reps_dic:
                        for entry.split(':')[0] in Reps_dic[key]:
                            
                            with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_Kef_sequences.faa", "a") as Kefout:
                                sequence  = seq_dic[entry.split(':')[0]][1]
                                id_tax =  seq_dic[entry.split(':')[0]][0]
                                Kef_list.append(id_tax.split('tax:')[1])
                                line = '>{}_Protein:Kef'.format(id_tax) + '\n' + str(sequence) + '\n'
                                Kefout.write(line)
                            Kefout.close()
                    else:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_Kef_sequences.faa", "a") as Kefout:
                            sequence  = seq_dic[key][1]
                            id_tax =  seq_dic[key][0]
                            Kef_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:Kef'.format(id_tax) + '\n' + str(sequence) + '\n'
                            Kefout.write(line)
                        Kefout.close()
            elif clade_dic[key] == '45068':
                    if key in Reps_dic:
                        for entry.split(':')[0] in Reps_dic[key]:
                            
                            with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_NhaA_sequences.faa", "a") as NhaAout:
                                sequence  = seq_dic[entry.split(':')[0]][1]
                                id_tax =  seq_dic[entry.split(':')[0]][0]
                                NhaA_list.append(id_tax.split('tax:')[1])
                                line = '>{}_Protein:NhaA'.format(id_tax) + '\n' + str(sequence) + '\n'
                                NhaAout.write(line)
                            NhaAout.close()
                    else:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_NhaA_sequences.faa", "a") as NhaAout:
                            sequence  = seq_dic[key][1]
                            id_tax =  seq_dic[key][0]
                            NhaA_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:NhaA'.format(id_tax) + '\n' + str(sequence) + '\n'
                            NhaAout.write(line)
                        NhaAout.close()
            elif clade_dic[key] == '57585':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_NhaS5_sequences.faa", "a") as NhaS5out:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            NhaS5_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:NhaS5'.format(id_tax) + '\n' + str(sequence) + '\n'
                            NhaS5out.write(line)
                        NhaS5out.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_NhaS5_sequences.faa", "a") as NhaS5out:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        NhaS5_list.append(id_tax.split('tax:')[1])

                        line = '>{}_Protein:NhaS5'.format(id_tax) + '\n' + str(sequence) + '\n'
                        NhaS5out.write(line)
                    NhaS5out.close()
            elif clade_dic[key] == '23026':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_Undescribed_CPA1_sequences.faa", "a") as Undescribed_CPA1out:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            Undescribed_CPA1_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:Undescribed_CPA1'.format(id_tax) + '\n' + str(sequence) + '\n'
                            Undescribed_CPA1out.write(line)
                        Undescribed_CPA1out.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_Undescribed_CPA1_sequences.faa", "a") as Undescribed_CPA1out:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        Undescribed_CPA1_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:Undescribed_CPA1'.format(id_tax) + '\n' + str(sequence) + '\n'
                        Undescribed_CPA1out.write(line)
                    Undescribed_CPA1out.close()
            elif clade_dic[key] == '11132':
                if key in Reps_dic:
                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_SOD2_sequences.faa", "a") as SOD2out:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            SOD2_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:SOD2'.format(id_tax) + '\n' + str(sequence) + '\n'
                            SOD2out.write(line)
                        SOD2out.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_SOD2_sequences.faa", "a") as  SOD2out:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        SOD2_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:SOD2'.format(id_tax) + '\n' + str(sequence) + '\n'
                        SOD2out.write(line)
                    SOD2out.close()
            elif clade_dic[key] == '14956':
                if key in Reps_dic:
                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_NhaPCPA1_sequences.faa", "a") as NhaPCPA1out:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            NhaP_CPA1_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:NhaPCPA1'.format(id_tax) + '\n' + str(sequence) + '\n'
                            NhaPCPA1out.write(line)
                        NhaPCPA1out.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_NhaPCPA1_sequences.faa", "a") as  NhaPCPA1out:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        NhaP_CPA1_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:NhaPCPA1'.format(id_tax) + '\n' + str(sequence) + '\n'
                        NhaPCPA1out.write(line)
                    NhaPCPA1out.close()
            elif clade_dic[key] == '8560':
                if key in Reps_dic:
                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_UndProkaryaCPA1IDK_sequences.faa", "a") as UndProkaryaCPA1IDKout:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            UndProkarya_CPA1_IDK_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:UndProkaryaCPA1IDK'.format(id_tax) + '\n' + str(sequence) + '\n'
                            UndProkaryaCPA1IDKout.write(line)
                        UndProkaryaCPA1IDKout.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_UndProkaryaCPA1IDK_sequences.faa", "a") as  UndProkaryaCPA1IDKout:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        UndProkarya_CPA1_IDK_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:UndProkaryaCPA1IDK'.format(id_tax) + '\n' + str(sequence) + '\n'
                        UndProkaryaCPA1IDKout.write(line)
                    UndProkaryaCPA1IDKout.close()
            elif clade_dic[key] == '53633':
                if key in Reps_dic:
                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_GerN_sequences.faa", "a") as GerNout:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            GerN_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:GerN'.format(id_tax) + '\n' + str(sequence) + '\n'
                            GerNout.write(line)
                        GerNout.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_GerN_sequences.faa", "a") as GerNout:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        GerN_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:GerN'.format(id_tax) + '\n' + str(sequence) + '\n'
                        GerNout.write(line)
                    GerNout.close()
            elif clade_dic[key] == '21649':
                if key in Reps_dic:
                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_CPA1_SL_sequences.faa", "a") as CPA1_SLout:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            CPA1_SL_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:CPA1_SL'.format(id_tax) + '\n' + str(sequence) + '\n'
                            CPA1_SLout.write(line)
                        CPA1_SLout.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_CPA1_SL_sequences.faa", "a") as CPA1_SLout:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]

                        CPA1_SL_list.append(id_tax.split('tax:')[1])
                        line = '>{}_Protein:CPA1_SL'.format(id_tax) + '\n' + str(sequence) + '\n'
                        CPA1_SLout.write(line)
                    CPA1_SLout.close()
            else:
                if key in Reps_dic:
                    for entry in Reps_dic[key]:
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_Uncharacterized_sequences.faa", "a") as Uncharacterizedout:
                            sequence  = seq_dic[entry.split(':')[0]][1]
                            id_tax =  seq_dic[entry.split(':')[0]][0]
                            Uncharacterized_list.append(id_tax.split('tax:')[1])
                            line = '>{}_Protein:Uncharacterized'.format(id_tax) + '\n' + str(sequence) + '\n'
                            Uncharacterizedout.write(line)
                        Uncharacterizedout.close()
                else:
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_april28/fl/Bac_Uncharacterized_sequences.faa", "a") as Uncharacterizedout:
                        sequence  = seq_dic[key][1]
                        id_tax =  seq_dic[key][0]
                        Uncharacterized_list.append(id_tax.split('tax:')[1])

                        line = '>{}_Protein:Uncharacterized'.format(id_tax) + '\n' + str(sequence) + '\n'
                        Uncharacterizedout.write(line)
                    Uncharacterizedout.close()

print(len(Kef_list))
print(len(CPA1_list))
print(len(NhaA_list))
print(len(Undescribed_CPA1_list))
print(len(UndProkarya_CPA1_IDK_list))
print(len(CHX_list))
print(len(GerN_list))
print(len(SOD2_list))
print(len(NhaS5_list))
print(len(Uncharacterized_list))
print(len(NhaP_CPA1_list))
print(len(CPA1_SL_list))

with open('/nesi/nobackup/uc04105/new_databases_May/final_april28/CPA_types_Bacteria_5may.tsv', 'w') as TAX_TYPE:
    header = 'GTDB_id' + '\t' + 'Kef' + '\t' + 'CPA1' + '\t' + 'NhaA' + '\t' + 'Undescribed_CPA1' + '\t' + 'NhaS5' + '\t' + 'CHX' + '\t' + 'SOD2' + '\t' + 'NhaP' + '\t' + 'CPA1_IDK' + '\t' + 'GerN' + '\t' + 'CPA1_SL' + '\t' + 'Uncharacaterized' + '\n'
    TAX_TYPE.write(header)
    total_count = 0
    for GTDB_tax in tax_list:
        line = GTDB_tax + '\t' + str(Kef_list.count(GTDB_tax)) + '\t'  + str(CPA1_list.count(GTDB_tax)) + '\t' + str(NhaA_list.count(GTDB_tax)) + '\t' + str(Undescribed_CPA1_list.count(GTDB_tax)) + '\t'  + str(NhaS5_list.count(GTDB_tax)) + '\t' + str(CHX_list.count(GTDB_tax)) + '\t'  + str(SOD2_list.count(GTDB_tax)) + '\t'  + str(NhaP_CPA1_list.count(GTDB_tax)) + '\t' + str(UndProkarya_CPA1_IDK_list.count(GTDB_tax)) + '\t' + str(GerN_list.count(GTDB_tax)) + '\t' + str(CPA1_SL_list.count(GTDB_tax)) + '\t' + str(Uncharacterized_list.count(GTDB_tax)) + '\n'
        TAX_TYPE.write(line)                   
