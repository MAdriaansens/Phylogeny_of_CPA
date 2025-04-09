import json
diction = '/nesi/nobackup/uc04105/Bacteria_hits_full_mar17.json'
with open(diction) as json_file:
        data = json.load(json_file)

final_list = []
intermediate = list(data.values())

for combi in intermediate:
    final_list.append(combi[1])
print(len(final_list))
 with open('/nesi/nobackup/uc04105/CPA_hitBacteria28032025.tsv', 'a') as B:
    header = 'accesion' + '\t'+ 'CheckM2_completeness' + '\t' + 'CheckM2_contamination' + '\t' + 'CPA_count' + '\t' + 'CPA_binary' + '\t'+ 'GTDB_taxonomy' + '\n'
    B.write(header)
    with open('/nesi/nobackup/uc04105/database/GTDB_220/bac120_metadata_r220.tsv', 'r') as M:
        for i in M:
            if (i.split('\t')[18]) == 't':
                accession = i.split('\t')[0]
                CheckM2_completeness = i.split('\t')[2]
                CheckM2_contamination = i.split('\t')[3]
                Taxonomy = i.split('\t')[19].replace(' ', '_')
                CPA_count = (final_list.count(Taxonomy))
                if CPA_count != 0:
                    CPA_binary = 1
                else:
                    CPA_binary = 0
                line = accession + '\t' + str(CheckM2_completeness) + '\t' + str(CheckM2_contamination) + '\t' + str(CPA_count) + '\t' + str(CPA_binary) + '\t' + Taxonomy + '\n'
                B.write(line)
