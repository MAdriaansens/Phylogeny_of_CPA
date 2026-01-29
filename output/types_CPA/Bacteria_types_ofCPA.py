
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
    if 'Bac' in key:
        if key[0].isnumeric():
            pass
        else:
            if clade_dic[key] == '35427':
                #now we unpack the representatives
                if key in Reps_dic:
                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        CPA1_dict[entry] = 'CPA1'
                        
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/CPA1_sequences.faa", "a") as CPA1out:
                            sequence  = seq_dic[entry]
                            line = '>{}_CPA1'.format(entry) + '\n' + str(sequence) + '\n'
                            CPA1out.write(line)
                        CPA1out.close()
                else:
                    CPA1_dict[key] = 'CPA1'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/CPA1_sequences.faa", "a") as CPA1out:
                        sequence  = seq_dic[key]
                        line = '>{}_CPA1'.format(key) + '\n' + str(sequence) + '\n'
                        CPA1out.write(line)
                    CPA1out.close()
    
            elif clade_dic[key] == '13164':
                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        Kef_dict[entry] = 'Kef'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Kef_sequences.faa", "a") as KEFout:
                            sequence  = seq_dic[entry]
                            line = '>{}_Kef'.format(entry) + '\n' + str(sequence) + '\n'
                            KEFout.write(line)
                        KEFout.close()
                else:
                    Kef_dict[key] = 'Kef'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Kef_sequences.faa", "a") as KEFout:
                        sequence  = seq_dic[key]
                        line = '>{}_Kef'.format(key) + '\n' + str(sequence) + '\n'
                        KEFout.write(line)
                    KEFout.close()
            elif clade_dic[key] == '4':
                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        
                        CPA2_dict[entry] = 'CPA2'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/CPA2_sequences.faa", "a") as CPA2out:
                            sequence  = seq_dic[entry]
                            line = '>{}_CPA2'.format(entry) + '\n' + str(sequence) + '\n'
                            CPA2out.write(line)
                        CPA2out.close()
                else:
                    CPA2_dict[key] = 'CPA2'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/CPA2_sequences.faa", "a") as CPA2out:
                        sequence  = seq_dic[key]
                        line = '>{}_CPA2'.format(key) + '\n' + str(sequence) + '\n'
                        CPA2out.write(line)
                    CPA2out.close()
            elif clade_dic[key] == '27962':
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        NhaA_dict[entry] = 'NhaA'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/NhaA_sequences.faa", "a") as NhaAout:
                            sequence  = seq_dic[entry]
                            line = '>{}_NhaA'.format(entry) + '\n' + str(sequence) + '\n'
                            NhaAout.write(line)
                        NhaAout.close()
                else:
                    NhaA_dict[key] = 'NhaA'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/NhaA_sequences.faa", "a") as NhaAout:
                        sequence  = seq_dic[key]
                        line = '>{}_NhaA'.format(key) + '\n' + str(sequence) + '\n'
                        NhaAout.write(line)
                    NhaAout.close()
            elif clade_dic[key] == '58708':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        Lates_dict[entry] = 'Latescibacteriota_cpa'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Latescibacteriota_cpa_sequences.faa", "a") as Latescibacteriota_out:
                            sequence  = seq_dic[entry]
                            line = '>{}_Latescibacteriota_cpa'.format(entry) + '\n' + str(sequence) + '\n'
                            Latescibacteriota_out.write(line)
                        Latescibacteriota_out.close()
                else:
                    Lates_dict[key] = 'Latescibacteriota_cpa'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Latescibacteriota_cpa_sequences.faa", "a") as Latescibacteriota_out:
                        sequence  = seq_dic[key]
                        line = '>{}_Latescibacteriota_cpa'.format(key) + '\n' + str(sequence) + '\n'
                        Latescibacteriota_out.write(line)
                    Latescibacteriota_out.close()

            elif clade_dic[key] == '58745':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        DxK_pseudo_dict[entry] = 'DxK_Pseudomonadota'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/DxK_Pseudomonadota_sequences.faa", "a") as DxK_Pseudomonadota_out:
                            sequence  = seq_dic[entry]
                            line = '>{}_DxK_Pseudomonadota'.format(entry) + '\n' + str(sequence) + '\n'
                            DxK_Pseudomonadota_out.write(line)
                        DxK_Pseudomonadota_out.close()
                else:
                    DxK_pseudo_dict[key] = 'DxK_Pseudomonadota'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/DxK_Pseudomonadota_sequences.faa", "a") as DxK_Pseudomonadota_out:
                        sequence  = seq_dic[key]
                        line = '>{}_DxK_Pseudomonadota'.format(key) + '\n' + str(sequence) + '\n'
                        DxK_Pseudomonadota_out.write(line)
                    DxK_Pseudomonadota_out.close()
            elif clade_dic[key] == '59314':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        NxKGamma_dict[entry] = 'NxK_Gammaproteobacteriota'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/NxK_Gammaproteobacteriota_sequences.faa", "a") as NxK_Gammaproteobacteriota_out:
                            sequence  = seq_dic[entry]
                            line = '>{}_NxK_Gammaproteobacteriota'.format(entry) + '\n' + str(sequence) + '\n'
                            NxK_Gammaproteobacteriota_out.write(line)
                        NxK_Gammaproteobacteriota_out.close()
                else:
                    NxKGamma_dict[key] = 'NxK_Gammaproteobacteriota'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/NxK_Gammaproteobacteriota_sequences.faa", "a") as NxK_Gammaproteobacteriota_out:
                        sequence  = seq_dic[key]
                        line = '>{}_NxK_Gammaproteobacteriota'.format(key) + '\n' + str(sequence) + '\n'
                        NxK_Gammaproteobacteriota_out.write(line)
                    NxK_Gammaproteobacteriota_out.close()
            elif clade_dic[key] == '65825':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        NhaS5_dict[entry] = 'NhS5'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/NhaS5_sequences.faa", "a") as NhaS5out:
                            sequence  = seq_dic[entry]
                            line = '>{}_NhaS5'.format(entry) + '\n' + str(sequence) + '\n'
                            NhaS5out.write(line)
                        NhaS5out.close()
                else:
                    NhaS5_dict[key] = 'NhaS5'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/NhaS5_sequences.faa", "a") as NhaS5out:
                        sequence  = seq_dic[key]
                        line = '>{}_NhaS5'.format(key) + '\n' + str(sequence) + '\n'
                        NhaS5out.write(line)
                    NhaS5out.close()
            elif clade_dic[key] == '65668':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        UncArc_dict[entry] = 'UncArc'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/UncARC_sequences.faa", "a") as UncARCout:
                            sequence  = seq_dic[entry]
                            line = '>{}_UncARC'.format(entry) + '\n' + str(sequence) + '\n'
                            UncARCout.write(line)
                        UncARCout.close()
                else:
                    UncArc_dict[key] = 'UncArc'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/UncARC_sequences.faa", "a") as UncARCout:
                        sequence  = seq_dic[key]
                        line = '>{}_UncARC'.format(key) + '\n' + str(sequence) + '\n'
                        UncARCout.write(line)
                    UncARCout.close()
            elif clade_dic[key] == '61608':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        UncProk_dict[entry] = 'Uncharacterized_Prokarya'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/UncPROK_sequences.faa", "a") as UncPROKout:
                            sequence  = seq_dic[entry]
                            line = '>{}_UncPROK'.format(entry) + '\n' + str(sequence) + '\n'
                            UncPROKout.write(line)
                        UncPROKout.close()
                else:
                    UncProk_dict[key] = 'Uncharcaterized_Prokarya'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/UncPROK_sequences.faa", "a") as UncPROKout:
                        sequence  = seq_dic[key]
                        line = '>{}_UncPROK'.format(key) + '\n' + str(sequence) + '\n'
                        UncPROKout.write(line)
                    UncPROKout.close()
            else:
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key].replace('[','').replace(']', '').replace('"', '').replace(' ', '').split(','):
                        Uncharacterized[entry] = 'Unc'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Unc_sequences.faa", "a") as Uncout:
                            sequence  = seq_dic[entry]
                            line = '>{}_Unc'.format(entry) + '\n' + str(sequence) + '\n'
                            Uncout.write(line)
                        Uncout.close()
                else:
                    Uncharacterized[key] = 'Unc'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Unc_sequences.faa", "a") as Uncout:
                        sequence  = seq_dic[key]
                        line = '>{}_Unc'.format(key) + '\n' + str(sequence) + '\n'
                        Uncout.write(line)
                    Uncout.close()
            
print(len(set(CPA1_dict.keys())))
print(len(set(CPA2_dict.keys())))  
print(len(set(Kef_dict.keys())))        
print(len(set(NhaA_dict.keys())))
print(len(set(Uncharacterized.keys())))
    
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
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Types_CPA_in_BacteriaGTDB226.tsv', 'w') as TAX_TYPE:
    header = 'GTDB_id' + '\t' + 'Uncharcterized' + '\t' + 'Kef' + '\t' + 'CPA1' + '\t' + 'CPA2' + '\t' + 'NhaA' + '\t' + 'NhaS5' + '\t' + 'Uncharacterized_Archaea' + '\t' + 'Uncharacterized_prokarya' + '\t' + 'DxK_Pseudomonadota' + '\t' + 'NxK_Gammaproteobacteriota' + '\t' + 'Latescibacteriota' + '\n'
    TAX_TYPE.write(header)
    total_count = 0
    for tax in Tax_list:
        GTDB_tax = tax
        line = GTDB_tax + '\t' + str(Unc_list_tax.count(GTDB_tax)) + '\t' + str(Kef_list_tax.count(GTDB_tax)) + '\t' + str(CPA1_list_tax.count(GTDB_tax)) + '\t' +  str(CPA2_list_tax.count(GTDB_tax)) + '\t' + str(NhaA_list_tax.count(GTDB_tax)) + '\t'+ str(NhaS5_list_tax.count(GTDB_tax))  + '\t' +  str(UncArc_list_tax.count(GTDB_tax)) + '\t' +str(UncProk_list_tax.count(GTDB_tax)) + '\t' +  str(DxK_pseudo_list_tax.count(GTDB_tax)) + '\t' +  str(NxKGamma_list_tax.count(GTDB_tax)) + '\t' + str(Lates_list_tax.count(GTDB_tax)) + '\n'
        TAX_TYPE.write(line)

