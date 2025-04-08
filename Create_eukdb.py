import os
from Bio import SeqIO


TSV_db = '/home/mad149/Eukarya_db.tsv'
Prot_dir = '/home/mad149/new_euk'
File_list = os.listdir(Prot_dir)
out_tsv = '/nesi/nobackup/uc04105/database/Euk_db_7April/Euk_db_7April_protein.tsv'
out_db = '/nesi/nobackup/uc04105/database/Euk_db_7April/Euk_db_7April.fasta'
len(File_list)

with open(out_tsv, 'a') as O:
    header = 'Protein_id' + '\t' + 'Species_name' + '\t' +  'original_prot_id' + '\t' + 'Larger_grouping' + '\t' + 'NCBI_taxonomy' + 'Sequence'
    O.write(header)
    with open(out_db, 'a') as F:
        prot_count = 0
        #internal count of proteins
        with open(TSV_db, 'r') as T:
            #count = 0
            next(T, None)
            for i in T:
                Download_link = i.split('\t')[-2]
                if '.gz' not in Download_link:
                    if 'protein.faa' in Download_link:
                        matching_file = Download_link.split('/')[-2] +'_protein.faa'
                    elif 'El_Paco' in Download_link:
                        matching_file = 'El_Paco_V3_proteins.fa'
                    elif 'okinawa_mozuku' in Download_link:
                        matching_file = 'okinawa_mozuku_S_ver2_protein_seq.fasta'
                    else:
                        matching_file = Download_link.split('/')[-1]
                else:
                    matching_file =Download_link.split('/')[-1].split('.gz')[0]
                Name = i.split('\t')[1] 
                for record in SeqIO.parse('{}/{}'.format(Prot_dir,matching_file), 'fasta'):
                    original_protein_id = record.id
                    sequence = str(record.seq)
                    Larger_grouping = i.split('\t')[-4]
                    Taxonomy = ';p_' + i.split('\t')[-5] + ';c__' + i.split('\t')[-6] + ';o_' + i.split('\t')[-7] + ';f_' + i.split('\t')[-8] + ';g_' + i.split('\t')[-9]
                    NCBI_ID = i.split('\t')[0]
                    Name = i.split('\t')[1]
                    Protein_id = 'EukA7_{}'.format(prot_count)
                    TLine = Protein_id + '\t' + Name + '\t' + original_protein_id + '\t' + Larger_grouping + '\t' + Taxonomy + '\t' + sequence + '\n'
                    O.write(TLine)
                    prot_count = prot_count + 1
                    FLine = '>' + Protein_id + '\n' + sequence + '\n' 
                    F.write(FLine)

                #print(Name, matching_file)
                #if matching_file in File_list:
                    #count = count + 1
                #else:
                    #print(matching_file)
        #print(count)

#test case

#test_list = []
#with open(out_tsv, 'r') as OT:
 #   next(OT, None)
  #  for Entry in OT:
   #     if Entry.split('\t')[1] not in test_list:
    #        test_list.append(Entry.split('\t')[1])
#print(test_list)
#print(len(test_list))
