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


Tax_dic = {}



Reps_dic = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Eukarya_clusters_7OKT.tsv', 'r') as Reps_Bac:
    next(Reps_Bac, None)
    for line in Reps_Bac:
        Rep = (line.split('\t')[0])
        Clusteroids = (line.split('\t')[-1].split('\n')[0])
        Reps_dic[Rep] = Clusteroids
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/fl_taxa/Eukarya/Eukarya_cpa_fl_taxa.json', 'r') as euk_data:
    Euk_dic = json.load(euk_data)

seq_dic = {}
Tax_dic = {}
Tax_list = []

for key in Euk_dic.keys():
    Tax_list.append(Euk_dic[key][0].split('_Taxa')[0])
    seq_dic[key] = Euk_dic[key][-1]
    Tax_dic[key] = Euk_dic[key][0]


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
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/CPA1_Euksequences.faa", "a") as CPA1out:
                             sequence  = seq_dic[entry]
                             line = '>{}_{}_CPA1'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                             CPA1out.write(line)
                        CPA1out.close()
                else:
                    CPA1_dict[key] = 'CPA1'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/CPA1_Euksequences.faa", "a") as CPA1out:
                             sequence  = seq_dic[key]
                             line = '>{}_{}_CPA1'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                             CPA1out.write(line)
                    CPA1out.close()
    
            elif clade_dic[key] == '13164':
                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        Kef_dict[entry] = 'Kef'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/Kef_Euksequences.faa", "a") as Kefout:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_Kef'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            Kefout.write(line)
                        Kefout.close()
                else:
                    Kef_dict[key] = 'Kef'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/Kef_Euksequences.faa", "a") as Kefout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_Kef'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        Kefout.write(line)
                    Kefout.close()
            elif clade_dic[key] == '4':
                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        CPA2_dict[entry] = 'CPA2'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/CPA2_Euksequences.faa", "a") as CPA2out:
                             sequence  = seq_dic[entry]
                             line = '>{}_{}_CPA2'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                             CPA2out.write(line)
                        CPA2out.close()
                else:
                    CPA2_dict[key] = 'CPA2'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/CPA2_Euksequences.faa", "a") as CPA2out:
                             sequence  = seq_dic[key]
                             line = '>{}_{}_CPA2'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                             CPA2out.write(line)
                    CPA2out.close()
            elif clade_dic[key] == '27962':
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        NhaA_dict[entry] = 'NhaA'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/NhaA_Euksequences.faa", "a") as Nhaout:
                             sequence  = seq_dic[entry]
                             line = '>{}_{}_NhaA'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                             Nhaout.write(line)
                        NhaAout.close()
                else:
                    NhaA_dict[key] = 'NhaA'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/NhaA_Euksequences.faa", "a") as Nhaout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_NhaA'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                        Nhaout.write(line)
                    Nhaout.close()
            elif clade_dic[key] == '58708':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        Lates_dict[entry] = 'Latescibacteriota_cpa'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/Latesci_Euksequences.faa", "a") as Latesciout:
                             sequence  = seq_dic[entry]
                             line = '>{}_{}_Latesci'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                             Latesciout.write(line)
                        Latesciout.close()
                else:
                    Lates_dict[key] = 'Latescibacteriota_cpa'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/Latesci_Euksequences.faa", "a") as Latesciout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_Latesci'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        Latesciout.write(line)
                    Latesciout.close()

            elif clade_dic[key] == '58745':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        DxK_pseudo_dict[entry] = 'DxK_Pseudomonadota'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/DxKPseudo_Euksequences.faa", "a") as DxK:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_DxKPseudo'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            DxK.write(line)
                        Dxk.close()
                else:
                    DxK_pseudo_dict[key] = 'DxK_Pseudomonadota'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/DxKPseudo_Euksequences.faa", "a") as DxK:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_DxKPseudo'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        DxK.write(line)
                    Dxk.close()
            elif clade_dic[key] == '59314':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        NxKGamma_dict[entry] = 'NxK_Gammaproteobacteriota'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/NxK_Euksequences.faa", "a") as NxK:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_NxK'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            NxK.write(line)
                        Nxk.close()
                else:
                    NxKGamma_dict[key] = 'NxK_Gammaproteobacteriota'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/NxK_Euksequences.faa", "a") as NxK:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_NxK'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        NxK.write(line)
                    Nxk.close()
            elif clade_dic[key] == '65825':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        NhaS5_dict[entry] = 'NhS5'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/NhaS5_Euksequences.faa", "a") as NhaS5:
                            sequence  = seq_dic[key]
                            line = '>{}_{}_NhaS5'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            NhaS5.write(line)
                        NhaS5.close()
                else:
                    NhaS5_dict[key] = 'NhaS5'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/NhaS5_Euksequences.faa", "a") as NhaS5:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_NhaS5'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        NhaS5.write(line)
                    NhaS5.close()
            elif clade_dic[key] == '65668':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        UncArc_dict[entry] = 'UncArc'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/UncArc_Euksequences.faa", "a") as UncArc:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_UncArc'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            UnArc.write(line)
                        UnArc.close()
                else:
                    UncArc_dict[key] = 'UncArc'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/UncArc_Euksequences.faa", "a") as UncArc:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_UncArc'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        UnArc.write(line)
                    UnArc.close()
            elif clade_dic[key] == '61608':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        UncProk_dict[entry] = 'Uncharacterized_Prokarya'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/UncPro_Euksequences.faa", "a") as UncPro:
                            sequence  = seq_dic[key]
                            line = '>{}_{}_UncPro'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            UncPro.write(line)
                        UncPro.close()
                else:
                    UncProk_dict[key] = 'Uncharcaterized_Prokarya'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/UncPro_Euksequences.faa", "a") as UncPro:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_UncPro'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        UncPro.write(line)
                    UncPro.close()
            else:
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace("[","").replace("]", "").replace("'", "").replace(" ", "").split(","):
                        Uncharacterized[entry] = 'Unc'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/Unc_Euksequences.faa", "a") as Unc:
                            sequence  = seq_dic[key]
                            line = '>{}_{}_Unc'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            Unc.write(line)
                        Unc.close()
                else:
                    Uncharacterized[key] = 'Unc'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Eukarya/Unc_Euksequences.faa", "a") as Unc:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_Unc'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        Unc.write(line)





            
print(len(set(CPA1_dict.keys())))
print(len(set(CPA2_dict.keys())))  
print(len(set(Kef_dict.keys())))        
print(len(set(NhaA_dict.keys())))
print(len(set(Uncharacterized.keys())))
    

#make sure this some is similar to grep -c '>' Archaea input filei

CPA1_list_tax = []
CPA1_tax_id = []
for key in CPA1_dict.keys():
    CPA1_list_tax.append(Tax_dic[key].split('_Taxa')[0])
    #WORK ON THIS AFTER NZMS
    #info_list = []
    #info_list.append(Tax_dic[key])
    #info_list.append('CPA1')



CPA2_list_tax = []
for key in CPA2_dict.keys():
    CPA2_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(CPA2_list_tax))

NhaA_list_tax = []
for key in NhaA_dict.keys():
    NhaA_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(NhaA_list_tax))

Kef_list_tax = []
for key in Kef_dict.keys():
    Kef_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(Kef_list_tax))


Unc_list_tax = []
for key in Uncharacterized.keys():
    Unc_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(Unc_list_tax))

Lates_list_tax = []
for key in Lates_dict.keys():
    Lates_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(Lates_list_tax))


DxK_pseudo_list_tax = []
for key in DxK_pseudo_dict.keys():
    DxK_pseudo_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(DxK_pseudo_list_tax))

NxKGamma_list_tax = []
for key in NxKGamma_dict.keys():
    NxKGamma_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(NxKGamma_list_tax))


NhaS5_list_tax = []
for key in NhaS5_dict.keys():
    NhaS5_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(NhaS5_list_tax))

UncArc_list_tax = []
for key in UncArc_dict.keys():
    UncArc_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(UncArc_list_tax))

UncProk_list_tax = []
for key in UncProk_dict.keys():
    UncProk_list_tax.append(Tax_dic[key].split('_Taxa')[0])
print(len(UncProk_list_tax))


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Types_CPA_in_Eukarya_GTDB226.tsv', 'w') as TAX_TYPE:
    header = 'species_name' + '\t' + 'Uncharcterized' + '\t' + 'Kef' + '\t' + 'CPA1' + '\t' + 'CPA2' + '\t' + 'NhaA' + '\t' + 'NhaS5' + '\t' + 'Uncharacterized_Archaea' + '\t' + 'Uncharacterized_prokarya' + '\t' + 'DxK_Pseudomonadota' + '\t' + 'NxK_Gammaproteobacteriota' + '\t' + 'Latescibacteriota' + '\n'
    TAX_TYPE.write(header)
    for tax in set(Tax_list):
        Organism_id = tax
        line = Organism_id+ '\t' + str(Unc_list_tax.count(tax)) + '\t' + str(Kef_list_tax.count(tax)) + '\t' + str(CPA1_list_tax.count(tax)) + '\t' +  str(CPA2_list_tax.count(tax)) + '\t' + str(NhaA_list_tax.count(tax)) + '\t'+ str(NhaS5_list_tax.count(tax)) + '\t' +  str(UncArc_list_tax.count(tax)) + '\t' +str(UncProk_list_tax.count(tax)) + '\t' +  str(DxK_pseudo_list_tax.count(tax)) + '\t' +  str(NxKGamma_list_tax.count(tax)) + '\t' + str(Lates_list_tax.count(tax)) + '\n'
        TAX_TYPE.write(line)

            
print(len(set(CPA1_dict.keys())))
print(len(set(CPA2_dict.keys())))  
print(len(set(Kef_dict.keys())))        
print(len(set(NhaA_dict.keys())))
print(len(set(Uncharacterized.keys())))    
