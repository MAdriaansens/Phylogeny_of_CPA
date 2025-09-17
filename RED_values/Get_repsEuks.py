Tax = []
same_species_count = 0
num_same_species_proteins = 0
count = 0
count_1=0
count_2=0
count_3=0
count_4=0
count_5=0
count_6=0
count_7=0
count_8=0
count_9=0
Large_list = ('Opisthokonta', 'Sar')
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Eukarya_clusters_16sept.tsv', 'w') as Out:
    Header = 'Representative_id' + '\t' + 'No_proteins_in_cluster' + '\t' + 'Major_taxa'  + '\t' + 'No_diff_species' +  '\t' + 'Full_taxonomy' + '\t' + 'Diff_tax' + '\t' + 'Shared_tax' + '\t' + 'Species_list' + '\t' + 'clusteroid_list' + '\n'
    Out.write(Header)
    print(Header)
    for euk in replist:
        Shared_tax = ''
        minor = []
        minor_1 = []
        minor_2 = []
        minor_3 = []
        minor_4 = []
        minor_5 = []
        minor_6 = []
        minor_7 = []
        minor_8 = []
        minor_9 = []
        minor_10 = []
        minor_11 = []
        minor_12 = []
        minor_13 = []
        minor_14 = []
        minor_15 = []
        Taxa_in_cluster_list = []
        clusteroid_list = Rep_dict[euk]
        Species_list = []
        Diff_tax = []
        Tax = []
        Major_taxa_list = []
        for prot in clusteroid_list:
            Taxa_in_cluster_list.append(id_tax_dic[prot].split(';p')[1])
            Species_list.append(id_tax_dic[prot].split('_;p')[0])
        if len(set(Taxa_in_cluster_list)) == 1:
            Shared_tax = 'same_species'
            Representative_id = euk
            No_proteins_in_cluster = (len(clusteroid_list))
            List_proteins = clusteroid_list
            Diff_tax = 'same species'
            No_diff_species = (len(set(Taxa_in_cluster_list)))
            Species_list = str(set(Species_list))
            Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[1])
            Full_taxonomy = Taxa_in_cluster_list[0].split(';c__')[1]
            if Major_taxa in Large_list:
                if 'Sar' == Major_taxa:
                    Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                elif 'Opisthokonta' == Major_taxa:
                    Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
            count += 1
            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + Major_taxa  + '\t' + str(No_diff_species) +  '\t' + Full_taxonomy + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
        else:

            #check if we have multiple diff major taxa
            for MJ in Taxa_in_cluster_list:
                Major_taxa_list.append(MJ.split(';c_')[1].split(';')[1])
            if len(set(Major_taxa_list)) > 1:
                print('multiple taxa')
                print(Major_taxa_list)
                raise TypeError('multiple taxa')            
            for i in Taxa_in_cluster_list:
                    minor.append(i.split(';')[4])
                    if i.split(';')[4] == 'Fungiincertaesedis':
                        minor_1.append('Fungi_incertae_sedis')
                    else:
                        minor_1.append(i.split(';')[4])

                    minor_2.append(i.split(';')[5])
                    minor_3.append(i.split(';')[6])
                    minor_4.append(i.split(';')[7])
                    if len(i.split(';')) > 8:
                        minor_5.append(i.split(';')[8])
                    if len(i.split(';')) > 9:
                        minor_6.append(i.split(';')[9])
                    if len(i.split(';')) > 10:
                        minor_7.append(i.split(';')[10])
                    if len(i.split(';')) > 11:
                        minor_8.append(i.split(';')[11])
                    if len(i.split(';')) > 12:
                        minor_9.append(i.split(';')[12])
            if len(set(minor)) != 1:

                Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[1])
                Representative_id = euk
                No_proteins_in_cluster = (len(clusteroid_list))
                List_proteins = clusteroid_list
                No_diff_species = (len(set(Taxa_in_cluster_list)))
                Species_list = str(set(Species_list))
                Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[1])
                Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] 
                Diff_tax = str(list(set(minor)))
                print(minor)
                if Major_taxa in Large_list:
                    if 'Sar' == Major_taxa:
                        Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                    elif 'Opisthokonta' == Major_taxa:
                        Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                count_1 += 1
                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
            else:
                if (len(set(minor_1))) != 1:
                    Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                    Representative_id = euk
                    No_proteins_in_cluster = (len(clusteroid_list))
                    List_proteins = clusteroid_list
                    No_diff_species = (len(set(Taxa_in_cluster_list)))
                    Species_list = str(set(Species_list))
                    Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                    Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] 
                    Diff_tax = str(list(set(minor_1)))
                    count_2 += 1
                    Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                else:
                    if (len(set(minor_2))) != 1:
                        Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[3])
                        Representative_id = euk
                        No_proteins_in_cluster = (len(clusteroid_list))
                        List_proteins = clusteroid_list
                        No_diff_species = (len(set(Taxa_in_cluster_list)))
                        Species_list = str(set(Species_list))
                        Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                        Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2]  + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3]
                        Diff_tax = str(list(set(minor_2)))
                        count_3 += 1
                        Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                    else:
                        if (len(set(minor_3))) != 1:
                            Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[4])
                            Representative_id = euk
                            No_proteins_in_cluster = (len(clusteroid_list))
                            List_proteins = clusteroid_list
                            No_diff_species = (len(set(Taxa_in_cluster_list)))
                            Species_list = str(set(Species_list))
                            Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                            Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[4]
                            Diff_tax = str(list(set(minor_3)))
                            count_4 += 1
                            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                        else:
                            if (len(set(minor_4))) != 1:
                                Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[5])
                                Representative_id = euk
                                No_proteins_in_cluster = (len(clusteroid_list))
                                List_proteins = clusteroid_list
                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                Species_list = str(set(Species_list))
                                Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                                Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[4] +  '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[5]
                                Diff_tax = str(list(set(minor_4)))
                                count_5 += 1
                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                            else:
                                if (len(set(minor_5))) != 1:
                                    Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[6])
                                    Representative_id = euk
                                    No_proteins_in_cluster = (len(clusteroid_list))
                                    List_proteins = clusteroid_list
                                    No_diff_species = (len(set(Taxa_in_cluster_list)))
                                    Species_list = str(set(Species_list))
                                    Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                                    Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[4] +  '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[5] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[6]
                                    Diff_tax = str(list(set(minor_5)))
                                    count_6 += 1
                                    Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                                else:
                                    if (len(set(minor_6))) != 1:
                                        Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[7])
                                        Representative_id = euk
                                        No_proteins_in_cluster = (len(clusteroid_list))
                                        List_proteins = clusteroid_list
                                        No_diff_species = (len(set(Taxa_in_cluster_list)))
                                        Species_list = str(set(Species_list))
                                        Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                                        Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[4] +  \
                                        '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[5] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[6] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[7]
                                        Diff_tax = str(list(set(minor_6)))
                                        count_7 += 1
                                        Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                                    else:
                                        if (len(set(minor_7))) != 1:
                                            Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[8])
                                            Representative_id = euk
                                            No_proteins_in_cluster = (len(clusteroid_list))
                                            List_proteins = clusteroid_list
                                            No_diff_species = (len(set(Taxa_in_cluster_list)))
                                            Species_list = str(set(Species_list))
                                            Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                                            Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[4] +  \
                                            '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[5] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[6] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[7] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[8]
                                            Diff_tax = str(list(set(minor_7)))
                                            count_8 += 1
                                            Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
                                           # print(Line)
                                        else:
                                            if (len(set(minor_8))) != 1:
                                                Shared_tax =str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[9])
                                                Representative_id = euk
                                                No_proteins_in_cluster = (len(clusteroid_list))
                                                List_proteins = clusteroid_list
                                                No_diff_species = (len(set(Taxa_in_cluster_list)))
                                                Species_list = str(set(Species_list))
                                                Major_taxa = str(Taxa_in_cluster_list[0].split(';c_')[1].split(';')[2])
                                                Full_taxonomy = 'Eukarya_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[1] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[2] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[3] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[4] +  \
                                                '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[5] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[6] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[7] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[8] + '_' + Taxa_in_cluster_list[0].split(';c__')[1].split(';')[9]
                                                Diff_tax = str(list(set(minor_8)))
                                                count_9 += 1
        
                                                Line = Representative_id + '\t' + str(No_proteins_in_cluster) + '\t' + str(Major_taxa)  + '\t' + str(No_diff_species) +  '\t' + str(Full_taxonomy) + '\t' + str(Diff_tax) + '\t' + Shared_tax + '\t' + Species_list + '\t' + str(clusteroid_list) + '\n'
        if len(Diff_tax) == 0:
            print(Line)
            print('potty')
        Out.write(Line)
print(count, count_1, count_2, count_3, count_4, count_5, count_6, count_7, count_8, count_9)
print(count+count_1+count_2+count_3 + count_4 + count_5 + count_6 + count_7 + count_8 + count_9)
