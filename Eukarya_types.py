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


Tax_dic = {}



Reps_dic = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Eukarya_clusters_7OKT.tsv', 'r') as Reps_Bac:
    next(Reps_Bac, None)
    for line in Reps_Bac:
        Rep = (line.split('\t')[0])
        Clusteroids = (line.split('\t')[-1].split('\n')[0])
        Reps_dic[Rep] = Clusteroids


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
    if 'Euk' in key:
        if key[0].isnumeric():
            pass
        else:
            if clade_dic[key] == '35427':
                #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        CPA1_dict[entry] = 'CPA1'
                else:
                    CPA1_dict[key] = 'CPA1'
                
    
            elif clade_dic[key] == '13164':
                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        Kef_dict[entry] = 'Kef'
                else:
                    Kef_dict[key] = 'Kef'
            elif clade_dic[key] == '4':
                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        CPA2_dict[entry] = 'CPA2'
                else:
                    CPA2_dict[key] = 'CPA2'
            elif clade_dic[key] == '27962':
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        NhaA_dict[entry] = 'NhaA'
                else:
                    NhaA_dict[key] = 'NhaA'
            elif clade_dic[key] == '58708':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        Lates_dict[entry] = 'Latescibacteriota_cpa'
                else:
                    Lates_dict[key] = 'Latescibacteriota_cpa'

            elif clade_dic[key] == '58745':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        DxK_pseudo_dict[entry] = 'DxK_Pseudomonadota'
                else:
                    DxK_pseudo_dict[key] = 'DxK_Pseudomonadota'
            elif clade_dic[key] == '59314':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        NxKGamma_dict[entry] = 'NxK_Gammaproteobacteriota'
                else:
                    NxKGamma_dict[key] = 'NxK_Gammaproteobacteriota'
            elif clade_dic[key] == '65825':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        NhaS5_dict[entry] = 'NhS5'
                else:
                    NhaS5_dict[key] = 'NhaS5'
            elif clade_dic[key] == '65668':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        UncArc_dict[entry] = 'UncArc'
                else:
                    UncArc_dict[key] = 'Uncharcaterized_Archaea'
            elif clade_dic[key] == '61608':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        UncProk_dict[entry] = 'Uncharacterized_Prokarya'
                else:
                    UncProk_dict[key] = 'Uncharcaterized_Prokarya'
            else:
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                       Uncharacterized[entry] = 'Unc'
                else:
                    Uncharacterized[key] = 'Unc'
            
print(len(set(CPA1_dict.keys())))
print(len(set(CPA2_dict.keys())))  
print(len(set(Kef_dict.keys())))        
print(len(set(NhaA_dict.keys())))
print(len(set(Uncharacterized.keys())))
    
from Bio import SeqIO
Tax_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/All_merged_PF00999.faa', 'fasta'):
    Tax_dic[record.id.split('_tax:')[0]] = record.id.split('_tax:')[1]
    Tax_list.append(record.id.split('_tax:')[1])

Tax_list = set(Tax_list)
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

