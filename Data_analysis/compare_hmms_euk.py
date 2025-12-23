euk_sequences_passed_thresholds = {}

from Bio import SeqIO

for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/Eukarya_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta','fasta'):
    euk_sequences_passed_thresholds[str(record.id)] = record.seq
Hmmsearch_Arc = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Eukarya/Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned_vsEuk.tsv'

Hmmsearch_pf00999 = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Eukarya/PF00999_vsEuk.tsv'


sequences_hit_by_pf009999 = []
with open(Hmmsearch_pf00999, 'r') as PF00999:
    for line in PF00999:
        if line[0] == '#':
            pass
        else:
            sequences_hit_by_pf009999.append(line.split(' ')[0])

            
sequences_hit_by_pf009999 = list(set(sequences_hit_by_pf009999))
print(len(sequences_hit_by_pf009999))


list_pf00999_passed = []
for sequence in sequences_hit_by_pf009999:
    if sequence in euk_sequences_passed_thresholds:
        list_pf00999_passed.append(sequence)
print(len(list_pf00999_passed))

sequences_hit_by_Arc = []
with open(Hmmsearch_Arc, 'r') as Arc:
    for line in Arc:
        if line[0] == '#':
            pass
        else:
            sequences_hit_by_Arc.append(line.split(' ')[0])

            
sequences_hit_by_Arc = list(set(sequences_hit_by_Arc))
print(len(sequences_hit_by_Arc))


list_Arc_passed = []
for sequence in sequences_hit_by_Arc:
    if sequence in euk_sequences_passed_thresholds:
        list_Arc_passed.append(sequence)
print(len(list_Arc_passed))


Hmmsearch_Bac = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Eukarya/Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned_vsEuk.tsv'


sequences_hit_by_Bac = []
with open(Hmmsearch_Bac, 'r') as Bac:
    for line in Bac:
        if line[0] == '#':
            pass
        else:
            sequences_hit_by_Bac.append(line.split(' ')[0])

            
sequences_hit_by_Bac = list(set(sequences_hit_by_Bac))
print(len(sequences_hit_by_Bac))


list_Bac_passed = []
for sequence in sequences_hit_by_Bac:
    if sequence in euk_sequences_passed_thresholds:
        list_Bac_passed.append(sequence)
print(len(list_Bac_passed))

Hmmsearch_Euk = '/nesi/nobackup/uc04105/new_databases_May/comparative_hmmsearch/Eukarya/Manual_seq_cov30_e05_seqid0.7_genafpair_alignedEukarya_vsEuk.tsv'


sequences_hit_by_Euk = []
with open(Hmmsearch_Euk, 'r') as Euk:
    for line in Euk:
        if line[0] == '#':
            pass
        else:
            sequences_hit_by_Euk.append(line.split(' ')[0])

            
sequences_hit_by_Euk = list(set(sequences_hit_by_Euk))
print(len(sequences_hit_by_Euk))


list_Euk_passed = []
for sequence in sequences_hit_by_Euk:
    if sequence in euk_sequences_passed_thresholds:
        list_Euk_passed.append(sequence)
print(len(list_Euk_passed))
