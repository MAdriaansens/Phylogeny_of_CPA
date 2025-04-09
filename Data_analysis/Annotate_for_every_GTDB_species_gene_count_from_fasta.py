import sys

#all sequences
all_sequences = sys.argv[1]
#in GTDB_protein.tsv
tsv = sys.argv[2]

from Bio import SeqIO
import time
start= time.time()

id_list = []
for record in SeqIO.parse(all_sequences, 'fasta'):
    id_list.append(record.id)

id_list.sort()
id_list = list(set(id_list))
length = len(id_list)

count = 0
dict_hits = {}
with open(tsv, 'r') as B:
    for i in B:
        name = i.split('\t')[0]
        if name in id_list:
            description = []
            description.append(i.split('\t')[-2])
            description.append(i.split('\t')[-3])
            dict_hits[name] = description
            count = count + 1
            numb_seconds = time.time() - start
            time_to_end = round(((numb_seconds/count)*length)-numb_seconds)
            print(count, length, ' ... Time passed: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))), 'Expected to finish in: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(time_to_end)))-1, int(time.strftime('%d', time.gmtime(time_to_end)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(time_to_end))), flush = True)  
B.close()
            
import json
with open("Bacteria_hits_full_Apr8_PF06450_full_run_MMseq.json", "w") as outfile:
    json.dump(dict_hits, outfile)
