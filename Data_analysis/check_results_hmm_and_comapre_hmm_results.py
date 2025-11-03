import os

hmmsearch_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/hmmsearch_for_scan/'


Arc_list =[]
Euk_list = []
Bac_list = []
PF00999_list = []
for output in os.listdir(hmmsearch_dir):
    if 'ARC' in output:
        with open('{}/{}'.format(hmmsearch_dir, output), 'r') as HMMscan:
            for line in HMMscan:
                if '#' in line:
                    pass
                else:

                    Arc_list.append(line.split(' ')[0])
    elif 'EUK' in output:
        with open('{}/{}'.format(hmmsearch_dir, output), 'r') as HMMscan:
            for line in HMMscan:
                if '#' in line:
                    pass
                else:
                    Euk_list.append(line.split(' ')[0])
    elif 'Bacteria' in output:
        with open('{}/{}'.format(hmmsearch_dir, output), 'r') as HMMscan:
            for line in HMMscan:
                if '#' in line:
                    pass
                else:
                    Bac_list.append(line.split(' ')[0])
    elif 'PF00999' in output:
        with open('{}/{}'.format(hmmsearch_dir, output), 'r') as HMMscan:
            for line in HMMscan:
                if '#' in line:
                    pass
                else:
                    PF00999_list.append(line.split(' ')[0])
print(len(Arc_list))
print(len(Euk_list))
print(len(Bac_list))
print(len(PF00999_list))
#the lists are all the sequences they have hit

#now we compare with the hmmalign, whom also has taxonomy
from Bio import SeqIO

hmmalign_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/cross_domain_final'

passed_hmmalign_dict = {}
for aligned_output in os.listdir(hmmalign_dir):
    if 'refilter' in aligned_output:
        for record in SeqIO.parse('{}/{}'.format(hmmalign_dir, aligned_output), 'fasta'):
            passed_hmmalign_dict[record.id.split('_tax:')[0]] = record.id.split('tax:')[1]
#full lengh hmmscan
subset_hit_domain = {}
hmmscan_fl_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/subset_sequence/Bacteria'
CPA_list = ('Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned','Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C')
for hmmscan_fl_output in os.listdir(hmmscan_fl_dir):
    if 'tsv' in hmmscan_fl_output:
            
        with open('{}/{}'.format(hmmscan_fl_dir, hmmscan_fl_output), 'r') as flscan:
            for line in flscan:
                if '#' in line:
                    pass
                else:
                    subset_hit_domain[line.split('_tax')[0].split(' ')[-1]] = line.split(' ')[0]


Arc_count = 0
Bac_count = 0
Euk_count = 0
PF00999_count = 0

Arc_hmmalign_list =[]
Euk_hmmalign_list = []
Bac_hmmalign_list = []
PF00999_hmmalign_list = []

for protein_id in Arc_list:
    if protein_id in passed_hmmalign_dict.keys():
        Arc_count += 1
        Arc_hmmalign_list.append(protein_id)
for protein_id in Euk_list:
    if protein_id in passed_hmmalign_dict:
        Euk_count += 1
        Euk_hmmalign_list.append(protein_id)

for protein_id in Bac_list:
    if protein_id in passed_hmmalign_dict:
        Bac_count += 1
        Bac_hmmalign_list.append(protein_id)

for protein_id in PF00999_list:
    if protein_id in passed_hmmalign_dict:
        PF00999_count += 1
        PF00999_hmmalign_list.append(protein_id)

#hte number of sequences hit whom pass the hmmalign step

print(Arc_count)
print(Bac_count)
print(Euk_count)
print(PF00999_count)   

domains_hit_Arc = []
domains_hit_Bac = []
domains_hit_Euk = []
domains_hit_PF00999 = []

for protein_id in Arc_hmmalign_list:
    if protein_id in subset_hit_domain.keys():
        domains_hit_Arc.append(subset_hit_domain[protein_id])

for protein_id in Bac_hmmalign_list:
    if protein_id in subset_hit_domain.keys():
        domains_hit_Bac.append(subset_hit_domain[protein_id])

for protein_id in Euk_hmmalign_list:
    if protein_id in subset_hit_domain.keys():
        domains_hit_Euk.append(subset_hit_domain[protein_id])
for protein_id in PF00999_hmmalign_list:
    if protein_id in subset_hit_domain.keys():
        domains_hit_PF00999.append(subset_hit_domain[protein_id])
from collections import Counter
Counter(domains_hit_Arc)



Counter(domains_hit_PF00999)
Counter(domains_hit_Bac)
Counter(domains_hit_Euk)

count = 0
count_x = 0
phyla_list = []
for protein_id in set(PF00999_hmmalign_list):
    if protein_id in subset_hit_domain.keys():
        if subset_hit_domain[protein_id] in CPA_list:
            count_x = count_x + 1
            if protein_id not in Euk_hmmalign_list:
                if protein_id not in Bac_hmmalign_list:
                    if protein_id not in Arc_hmmalign_list:
                        phyla_list.append(passed_hmmalign_dict[protein_id])
                        count =  count +1
print(count)
print(count_x)
#switching ou whcih hmmalign_lists are included allows for us to identify unique sequences


class_list = []
sub_phyla_list = []
order_list = []
family_list =[]
for phyla in phyla_list:
    sub_phyla_list.append(phyla.split(';p__')[1].split(';c')[0])
    class_list.append(phyla.split(';c__')[1].split(';o')[0])
    order_list.append(phyla.split(';o__')[1].split(';f')[0])
    family_list.append(phyla.split(';f__')[1].split(';g')[0])

print(Counter(order_list))
