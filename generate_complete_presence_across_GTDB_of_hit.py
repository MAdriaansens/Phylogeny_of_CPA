metadata = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_metadata_r220.tsv'
list_completeness = []
count = 0
miss_tax_list = []
with open('complete.tsv', 'a') as C:
    header = 'accession' + '\t' + 'completeness' + '\n'
    #C.write(header)
    with open(metadata, 'r') as Meta:
        for line in Meta:
            if line.split('\t')[18] != 't':
                pass
            else:
                if (line.split('\t')[19].replace(' ', '_')) not in Uniq_tax_list:
                    if float(line.split('\t')[2]) > 80:
                        list_completeness.append(line.split('\t')[19].split(';o__')[1].split(';f')[0])
Counter(list_completeness)
 metadata = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_metadata_r220.tsv'
list_completeness = []
count = 0
miss_tax_list = []
with open('complete.tsv', 'a') as C:
    header = 'accession' + '\t' + 'completeness' + '\n'
    #C.write(header)
    with open(metadata, 'r') as Meta:
        for line in Meta:
            if line.split('\t')[18] != 't':
                pass
            else:
                if (line.split('\t')[19].replace(' ', '_')) not in Uniq_tax_list:
                    if float(line.split('\t')[2]) > 80:
                        list_completeness.append(line.split('\t')[19].split(';o__')[1].split(';f')[0])
Counter(list_completeness)

metadata = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_metadata_r220.tsv'
full_cluster = '/nesi/project/uc04105/Archaea_hits_full_mar14.json'
hit_tax = []
import json

with open(full_cluster) as json_file:
    data = json.load(json_file)

for i in list(data.values()):
    hit_tax.append(i[1])
print(len(hit_tax))


tax_list =[]
for i in data:
    tax_list.append(data[i][1])
    Uniq_tax_list = list(tax_list)



total_list = []
for i in tax_list:
    print(i)
    break
Phyla_tax_list = []
Class_tax_list = []
Order_tax_list = []
Family_tax_list = []
Genera_tax_list = []
Species_tax_list = []
for entry in tax_list:
    Phyla_tax_list.append(entry.split(';p__')[1].split(';c')[0])
    Class_tax_list.append(entry.split(';c__')[1].split(';o')[0])

    Order_tax_list.append(entry.split(';o__')[1].split(';f')[0])

    Family_tax_list.append(entry.split(';f__')[1].split(';g')[0])
    Genera_tax_list.append(entry.split(';g__')[1].split(';s')[0])
    Species_tax_list.append(entry.split(';s__')[1])
    if 'Asgardarchaeota' == entry.split(';p__')[1].split(';c')[0]:
        total_list.append(entry)



from collections import Counter
Counter(total_list)

