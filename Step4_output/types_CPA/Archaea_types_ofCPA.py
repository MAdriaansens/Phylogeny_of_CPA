

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
                        
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcCPA1_sequences.faa", "a") as CPA1out:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:CPA1'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            CPA1out.write(line)
                        CPA1out.close()
                else:
                    CPA1_dict[key] = 'CPA1'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcCPA1_sequences.faa", "a") as CPA1out:
                        sequence  = seq_dic[key]
    
                        line = '>{}_{}_protein:CPA1'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        CPA1out.write(line)
                    CPA1out.close()
    
            elif clade_dic[key] == '13164':
                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        Kef_dict[entry] = 'Kef'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcKef_sequences.faa", "a") as KEFout:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:Kef'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            KEFout.write(line)
                        KEFout.close()
                else:
                    Kef_dict[key] = 'Kef'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcKef_sequences.faa", "a") as KEFout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:Kef'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        KEFout.write(line)
                    KEFout.close()
            elif clade_dic[key] == '4':
                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        
                        CPA2_dict[entry] = 'CPA2'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcCPA2_sequences.faa", "a") as CPA2out:
                            sequence  = seq_dic[entry]

                            line = '>{}_{}_protein:CPA2'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            CPA2out.write(line)
                        CPA2out.close()
                else:
                    CPA2_dict[key] = 'CPA2'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcCPA2_sequences.faa", "a") as CPA2out:
                        sequence  = seq_dic[key]

                        line = '>{}_{}_protein:CPA2'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        CPA2out.write(line)
                    CPA2out.close()
            elif clade_dic[key] == '27962':
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        NhaA_dict[entry] = 'NhaA'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcNhaA_sequences.faa", "a") as NhaAout:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:NhaA'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            NhaAout.write(line)
                        NhaAout.close()
                else:
                    NhaA_dict[key] = 'NhaA'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcNhaA_sequences.faa", "a") as NhaAout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:NhaA'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        NhaAout.write(line)
                    NhaAout.close()
            elif clade_dic[key] == '58708':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        Lates_dict[entry] = 'Latescibacteriota_cpa'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcLatescibacteriota_cpa_sequences.faa", "a") as Latescibacteriota_out:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:Latescibacteriota_cpa'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            Latescibacteriota_out.write(line)
                        Latescibacteriota_out.close()
                else:
                    Lates_dict[key] = 'Latescibacteriota_cpa'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcLatescibacteriota_cpa_sequences.faa", "a") as Latescibacteriota_out:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:Latescibacteriota_cpa'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        Latescibacteriota_out.write(line)
                    Latescibacteriota_out.close()

            elif clade_dic[key] == '58745':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        DxK_pseudo_dict[entry] = 'DxK_Pseudomonadota'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcDxK_Pseudomonadota_sequences.faa", "a") as DxK_Pseudomonadota_out:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:DxK_Pseudomonadota'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            DxK_Pseudomonadota_out.write(line)
                        DxK_Pseudomonadota_out.close()
                else:
                    DxK_pseudo_dict[key] = 'DxK_Pseudomonadota'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcDxK_Pseudomonadota_sequences.faa", "a") as DxK_Pseudomonadota_out:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:DxK_Pseudomonadota'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        DxK_Pseudomonadota_out.write(line)
                    DxK_Pseudomonadota_out.close()
            elif clade_dic[key] == '59314':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        NxKGamma_dict[entry] = 'NxK_Gammaproteobacteriota'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcNxK_Gammaproteobacteriota_sequences.faa", "a") as NxK_Gammaproteobacteriota_out:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:NxK_Gammaproteobacteriota'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            NxK_Gammaproteobacteriota_out.write(line)
                        NxK_Gammaproteobacteriota_out.close()
                else:
                    NxKGamma_dict[key] = 'NxK_Gammaproteobacteriota'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcNxK_Gammaproteobacteriota_sequences.faa", "a") as NxK_Gammaproteobacteriota_out:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:NxK_Gammaproteobacteriota'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        NxK_Gammaproteobacteriota_out.write(line)
                    NxK_Gammaproteobacteriota_out.close()
            elif clade_dic[key] == '65825':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        NhaS5_dict[entry] = 'NhS5'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcNhaS5_sequences.faa", "a") as NhaS5out:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:NhaS5'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            NhaS5out.write(line)
                        NhaS5out.close()
                else:
                    NhaS5_dict[key] = 'NhaS5'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/ArcNhaS5_sequences.faa", "a") as NhaS5out:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:NhaS5'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        NhaS5out.write(line)
                    NhaS5out.close()
            elif clade_dic[key] == '65668':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        UncArc_dict[entry] = 'UncArc'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Arc_UncARC_sequences.faa", "a") as UncARCout:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:UncARC'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            UncARCout.write(line)
                        UncARCout.close()
                else:
                    UncArc_dict[key] = 'UncArc'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Arc_UncARC_sequences.faa", "a") as UncARCout:
                        sequence  = seq_dic[key]
                        line = '>{}_protein:UncARC'.format(key, Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        UncARCout.write(line)
                    UncARCout.close()
            elif clade_dic[key] == '61608':
                                                            #now we unpack the representatives
                if key in Reps_dic:

                    for entry in Reps_dic[key]:
                        UncProk_dict[entry] = 'Uncharacterized_Prokarya'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Arc_UncPROK_sequences.faa", "a") as UncPROKout:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:UncPROK'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            UncPROKout.write(line)
                        UncPROKout.close()
                else:
                    UncProk_dict[key] = 'Uncharcaterized_Prokarya'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Arc_UncPROK_sequences.faa", "a") as UncPROKout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:UncPROK'.format(key,Tax_dic[key]) + '\n' + str(sequence) + '\n'
                        UncPROKout.write(line)
                    UncPROKout.close()
            else:
                                                            #now we unpack the representatives
                if key in Reps_dic:
                    
                    for entry in Reps_dic[key]:
                        Uncharacterized[entry] = 'Unc'
                        with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Arc_Unc_sequences.faa", "a") as Uncout:
                            sequence  = seq_dic[entry]
                            line = '>{}_{}_protein:Unc'.format(entry, Tax_dic[entry]) + '\n' + str(sequence) + '\n'
                            Uncout.write(line)
                        Uncout.close()
                else:
                    Uncharacterized[key] = 'Unc'
                    with open("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/fl_sequences_types/Arc_Unc_sequences.faa", "a") as Uncout:
                        sequence  = seq_dic[key]
                        line = '>{}_{}_protein:Unc'.format(key,Tax_dic[key]) + '\n' + str(sequence) + '\n'
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
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Types_CPA_in_ArchaeaGTDB226.tsv', 'w') as TAX_TYPE:
    header = 'GTDB_id' + '\t' + 'Uncharcterized' + '\t' + 'Kef' + '\t' + 'CPA1' + '\t' + 'CPA2' + '\t' + 'NhaA' + '\t' + 'NhaS5' + '\t' + 'Uncharacterized_Archaea' + '\t' + 'Uncharacterized_prokarya' + '\t' + 'DxK_Pseudomonadota' + '\t' + 'NxK_Gammaproteobacteriota' + '\t' + 'Latescibacteriota' + '\n'
    TAX_TYPE.write(header)
    total_count = 0
    for tax in Tax_list:
        GTDB_tax = tax
        line = GTDB_tax + '\t' + str(Unc_list_tax.count(GTDB_tax)) + '\t' + str(Kef_list_tax.count(GTDB_tax)) + '\t' + str(CPA1_list_tax.count(GTDB_tax)) + '\t' +  str(CPA2_list_tax.count(GTDB_tax)) + '\t' + str(NhaA_list_tax.count(GTDB_tax)) + '\t'+ str(NhaS5_list_tax.count(GTDB_tax))  + '\t' +  str(UncArc_list_tax.count(GTDB_tax)) + '\t' +str(UncProk_list_tax.count(GTDB_tax)) + '\t' +  str(DxK_pseudo_list_tax.count(GTDB_tax)) + '\t' +  str(NxKGamma_list_tax.count(GTDB_tax)) + '\t' + str(Lates_list_tax.count(GTDB_tax)) + '\n'
        TAX_TYPE.write(line)
