new_json = '/nesi/nobackup/uc04105/Bacteria_hits_full_Apr8_PF03600_BlankHMM_MMseq.json'
old_tsv = '/nesi/nobackup/uc04105/CPA_hitBacteria28032025.tsv'
new_tsv = '/nesi/nobackup/uc04105/CPA_hitBacteria28032025_NhaD_edit.tsv'

import json

with open(new_json) as json_file:
    data = json.load(json_file)
    
hit_tax = []
for i in list(data.values()):
    hit_tax.append(i[0])

with open(new_tsv, 'a') as N:
    with open(old_tsv, 'r') as O:
        for entry in O:
            header =entry.split('\n')[0] + '\t' + 'NhaD_count' +'\t'+ 'NhaD_binary'+  '\n'
            N.write(header)
            break
        next(O, None)
        for row in O:
            NhaD_count = (hit_tax.count(row.split('\t')[0]))

            if NhaD_count != 0:
                edit_row =  row.split('\n')[0] + '\t' + str(NhaD_count) + '\t' + '1' + '\n'
                N.write(edit_row)
            else:
                edit_row =  row.split('\n')[0] + '\t' + '0' + '\t' + '0' + '\n'
                N.write(edit_row)
