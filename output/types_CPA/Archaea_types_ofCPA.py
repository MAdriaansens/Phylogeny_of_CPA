import json

clade_dic = {}

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Red_informed_clade_20_lg_cat_gamma_RED_interval.tsv', 'r') as LG_Cat_Gamma_clade:
    next(LG_Cat_Gamma_clade, None)
    for line in LG_Cat_Gamma_clade:
        protein_id = line.split('\t')[0]
        domain = line.split('\t')[1]
        clade_dic[protein_id] = domain
domain_list = []
for key in clade_dic.keys():
    domain_list.append(clade_dic[key])

from collections import Counter
print(Counter(domain_list))

#this matches the family to which a termini belongs
Tax_dic = {}
Tax_list = []

import os

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/MMseq/All_Arc_CPAseq_GTDB226_ids_tax.json')as json_file:
    d = json.load(json_file)
    for key in d.keys():
        Protein_id =  key

        Taxonomy  = d[key].split('__GTDB_tax:')[0].split('_id')[1]
        Tax_dic[Protein_id] = Taxonomy
        Tax_list.append(Taxonomy)
Tax_list = set(Tax_list)

#print uniq taxa with CPA, should be 5551 species. 
print(len(Tax_list))

#print total number of sequences, should be equal to 14101 (total number of CPAs)
print(len(Tax_dic.keys()))

from Bio import SeqIO

All_ids = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Archaea_passed_all_filters_7OKT_GTDBreps_alignedPF00999.fasta', 'fasta'):
    All_ids.append(record.id)
Reps_dic = {}
Clustered_sequences_list = []
#takes all clusteroids
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_values_Archaea_clusters_7OKT.tsv', 'r') as Reps_Arc:
    total_sequences = 0
    next(Reps_Arc, None)
    for line in Reps_Arc:
        Rep = (line.split('\t')[0])
        Clusteroids = list(line.split('\t')[-3].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","))
        Clustered_sequences_list.extend(Clusteroids)
        Reps_dic[Rep] = Clusteroids
Total_seq= 0


for entry in All_ids: 
    if entry not in Clustered_sequences_list:
        inlist = []
        inlist.append(entry)
        Reps_dic[entry] = inlist

print(len(Reps_dic.keys()))
#print(len(set(Reps_dic.keys())))


for key in Reps_dic.keys():
    Total_seq = Total_seq +len(Reps_dic[key])
print(Total_seq)







CPA1_dict = {}
CPA2_dict = {}
Kef_dict = {}
NhaA_dict = {}
Lates_dict={}
Uncharacterized  = {}
DxK_pseudo_dict = {}
NxKGamma_dict = {}
NhaS5_dict = {}
UncArc_dict = {}
UncProk_dict = {}


for key in clade_dic.keys():
    if 'Arc' in key:
        if key[0].isnumeric():
            pass
        else:
            if clade_dic[key] == '35427':
                #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        CPA1_dict[entry] = 'CPA1'
                else:
                    CPA1_dict[key] = 'CPA1'
                
    
            elif clade_dic[key] == '13164':
                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        Kef_dict[entry] = 'Kef'
                else:
                    Kef_dict[key] = 'Kef'
            elif clade_dic[key] == '4':
                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        CPA2_dict[entry] = 'CPA2'
                else:
                    CPA2_dict[key] = 'CPA2'
            elif clade_dic[key] == '27962':
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        NhaA_dict[entry] = 'NhaA'
                else:
                    NhaA_dict[key] = 'NhaA'
            elif clade_dic[key] == '58708':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        Lates_dict[entry] = 'Latescibacteriota_cpa'
                else:
                    Lates_dict[key] = 'Latescibacteriota_cpa'

            elif clade_dic[key] == '58745':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        DxK_pseudo_dict[entry] = 'DxK_Pseudomonadota'
                else:
                    DxK_pseudo_dict[key] = 'DxK_Pseudomonadota'
            elif clade_dic[key] == '59314':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        NxKGamma_dict[entry] = 'NxK_Gammaproteobacteriota'
                else:
                    NxKGamma_dict[key] = 'NxK_Gammaproteobacteriota'
            elif clade_dic[key] == '65825':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        NhaS5_dict[entry] = 'NhS5'
                else:
                    NhaS5_dict[key] = 'NhaS5'
            elif clade_dic[key] == '65668':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        UncArc_dict[entry] = 'UncArc'
                else:
                    UncArc_dict[key] = 'Uncharcaterized_Archaea'
            elif clade_dic[key] == '61608':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        UncProk_dict[entry] = 'Uncharacterized_Prokarya'
                else:
                    UncProk_dict[key] = 'Uncharcaterized_Prokarya'
            else:
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                       Uncharacterized[entry] = 'Unc'
                else:
                    Uncharacterized[key] = 'Unc'
            
print(len(set(CPA1_dict.keys())) + (len(set(CPA2_dict.keys()))) + (len(set(Kef_dict.keys()))) + (len(set(NhaA_dict.keys()))) + (len(set(Uncharacterized.keys()))) + (len(set(NxKGamma_dict.keys()))) + (len(set(NhaS5_dict.keys()))) + (len(set(UncArc_dict.keys()))) + (len(set(UncProk_dict.keys())))+ (len(set(DxK_pseudo_dict.keys()))))

#make sure this some is similar to grep -c '>' Archaea input filei
CPA1_list_tax = []
CPA1_tax_id = []
for key in CPA1_dict.keys():
    CPA1_list_tax.append(Tax_dic[key])
    #WORK ON THIS AFTER NZMS
    #info_list = []
    #info_list.append(Tax_dic[key])
    #info_list.append('CPA1')
print(len(CPA1_list_tax))
print(CPA1_list_tax[1])

CPA2_list_tax = []
for key in CPA2_dict.keys():
    CPA2_list_tax.append(Tax_dic[key])
print(len(CPA2_list_tax))

NhaA_list_tax = []
for key in NhaA_dict.keys():
    NhaA_list_tax.append(Tax_dic[key])
print(len(NhaA_list_tax))

Kef_list_tax = []
for key in Kef_dict.keys():
    Kef_list_tax.append(Tax_dic[key])
print(len(Kef_list_tax))


Unc_list_tax = []
for key in Uncharacterized.keys():
    Unc_list_tax.append(Tax_dic[key])
print(len(Unc_list_tax))

Lates_list_tax = []
for key in Lates_dict.keys():
    Lates_list_tax.append(Tax_dic[key])
print(len(Lates_list_tax))


DxK_pseudo_list_tax = []
for key in DxK_pseudo_dict.keys():
    DxK_pseudo_list_tax.append(Tax_dic[key])
print(len(DxK_pseudo_list_tax))

NxKGamma_list_tax = []
for key in NxKGamma_dict.keys():
    NxKGamma_list_tax.append(Tax_dic[key])
print(len(NxKGamma_list_tax))


NhaS5_list_tax = []
for key in NhaS5_dict.keys():
    NhaS5_list_tax.append(Tax_dic[key])
print(len(NhaS5_list_tax))

UncArc_list_tax = []
for key in UncArc_dict.keys():
    UncArc_list_tax.append(Tax_dic[key])
print(len(UncArc_list_tax))

UncProk_list_tax = []
for key in UncProk_dict.keys():
    UncProk_list_tax.append(Tax_dic[key])
print(len(UncProk_list_tax))


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Types_CPA_in_ArchaeaGTDB226.tsv', 'w') as TAX_TYPE:
    header = 'GTDB_id' + '\t' + 'Uncharcterized' + '\t' + 'Kef' + '\t' + 'CPA1' + '\t' + 'CPA2' + '\t' + 'NhaA' + '\t' + 'NhaS5' + '\t' + 'Uncharacterized_Archaea' + '\t' + 'Uncharacterized_prokarya' + '\t' + 'DxK_Pseudomonadota' + '\t' + 'NxK_Gammaproteobacteriota' + '\t' + 'Latescibacteriota' + '\n'
    TAX_TYPE.write(header)
    total_count = 0
    for tax in Tax_list:
        GTDB_tax = tax
        line = GTDB_tax + '\t' + str(Unc_list_tax.count(GTDB_tax)) + '\t' + str(Kef_list_tax.count(GTDB_tax)) + '\t' + str(CPA1_list_tax.count(GTDB_tax)) + '\t' +  str(CPA2_list_tax.count(GTDB_tax)) + '\t' + str(NhaA_list_tax.count(GTDB_tax)) + '\t'+ str(NhaS5_list_tax.count(GTDB_tax)) + '\t' +  str(UncArc_list_tax.count(GTDB_tax)) + '\t' +str(UncProk_list_tax.count(GTDB_tax)) + '\t' +  str(DxK_pseudo_list_tax.count(GTDB_tax)) + '\t' +  str(NxKGamma_list_tax.count(GTDB_tax)) + '\t' + str(Lates_list_tax.count(GTDB_tax)) + '\n'
        TAX_TYPE.write(line)
