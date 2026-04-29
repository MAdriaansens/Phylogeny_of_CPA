import sys
passed_seq = sys.argv[1] #'/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/results/HMMalign/PF00999/Eukfoundseq_vsEukarya_subset1_e03_mmseq_alignedPF00999.fasta.fasta'
tsv_file  = sys.argv[2] #'/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/Euk_DB/tsv/Euk_db_May_protein_chunk_1.tsv'
output = sys.argv[3] #'/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/results/HMMalign/PF00999/Eukfoundseq_vsEukarya_subset1_e03_mmseq_alignedPF00999_FL_nogrep.fasta'


from Bio import SeqIO

seq_list = []

for record in SeqIO.parse(passed_seq, 'fasta'):
    seq_list.append(record.id.split('_tax')[0])
unique_list = list(set(seq_list))
print(unique_list[1])


with open("{}".format(tsv_file)) as infile:
    file1 = open("{}".format(output), "w")  # append mode
    count = 0
    for line in infile:
        entry = line.split("\t")
        if entry[0] not in unique_list:
            pass
        else:
            if len(entry) < 5:
                pass
            else:
                sequence = ('>'+ entry[0] + '_tax:' +  entry[1] + '\n' + entry[-1])
                file1.write(sequence)
                count = count + 1

    file1.close()
