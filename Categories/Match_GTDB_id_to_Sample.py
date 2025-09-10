Cat_dict={}
#Load the manually edited file
with open('/nesi/nobackup/uc04105/new_databases_May/Sample_sites - Bacteria categories.tsv', 'r') as Categories:
    next(Categories, None)
    for Category in Categories:
        Sample = Category.split('\t')[0]
        Cat = Category.split('\t')[-1].split('\n')[0]
        Cat_dict[Sample]=Cat

  with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/IT_Bacteria_9SEPT.tsv', 'r') as R:
    next(R, None)
    with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_GTDB226_categories.tsv', 'w') as Ct:
        header = 'GTDB_id' + '\t' + 'Category' + '\n'
        Ct.write(header)
        for line in R:
            GTDB_id = (line.split('\t')[0])
            Sample = line.split('\t')[4].lower()
            if Sample != 'none':
                Category = Cat_dict[Sample]
                print(Category)
            else:
                Category = 'none described'
            Wline = GTDB_id + '\t' + Category + '\n'
            Ct.write(Wline)
