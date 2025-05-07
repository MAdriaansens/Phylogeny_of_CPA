foi = 'PF00999_07_vsEuk.tsv'
Value_dict = {}
with open(foi, 'r') as F:
    for line in F:
        i = line.split(" ")
        if '#' not in i[0]:
            list_value = []
            entry = i[0]
            bitscore =(i[30])
            for value in i:
                if '.' in value:
                    list_value.append(value)
            Value_dict[entry] = list_value

e_value_distribution = 'PF00999_07_vsEuk_distribution_evalue.tsv'
with open(e_value_distribution, 'a') as E:
    E.write('protein_id' +  '\t' + 'E-value' + '\t' + 'bitscore' + '\n')
    for key in Value_dict.keys():
        value = Value_dict[key]
        if len(value) == 8:
            #normal values
            evalue = value[1]
            bitscore = value[2]
        elif len(value) == 6:
            evalue = '0.0011'
            #note that it is greater most likely
            bitscore = value[1]
        elif len(value) == 7:
            if 'e' in value[1]:
                evalue = value[1]
                bitscore = value[2]
            elif 'e' in value[3]:
                evalue = value[3]
                bitscore = value[1]
        line = key + '\t' + evalue + '\t' + bitscore + '\n'
        E.write(line)
