import os

Arc_dir = '/nesi/nobackup/uc04105/Likelihood_CPA/Archaea'
Bac_dir = '/nesi/nobackup/uc04105/Likelihood_CPA/Bacteria'
Euk_dir = '/nesi/nobackup/uc04105/Likelihood_CPA/Eukarya'
Seed_dir = '/nesi/nobackup/uc04105/Likelihood_CPA/Seeds'

AIK_Arc = '{}/AIK_Arc.txt'.format(Arc_dir)
BIC_Arc = '{}/BIC_Arc.txt'.format(Arc_dir)

AIC_dic = {}
BIC_dic = {}

#we make a dic matching the model chosen for each number of subset of sequences 

with open(AIK_Arc, 'r') as AIK:
    for line in AIK:
        model_chosen = line.split(': ')[-1].split(' ')[0]
        number_of_file = line.split('.fasta.log')[0].split('_')[-1]
        AIC_dic[number_of_file] = model_chosen
        

with open(BIC_Arc, 'r') as BIC:
    for line in BIC:
        model_chosen = line.split(': ')[-1].split(' ')[0]
        number_of_file = line.split('.fasta.log')[0].split('_')[-1]
        BIC_dic[number_of_file] = model_chosen

count = 0

with open('/nesi/nobackup/uc04105/Likelihood_CPA/Archaea/Arc_likelihood_CPA_random100.tsv', 'w') as random100:
    header = 'random_set_number' + '\t'  + 'AIC_model' + '\t' + 'AIC_score' + '\t' + 'df_A' + '\t' +  + 'BIC_model' + '\t' + 'BIC_score' + '\t' + 'df_B' + '\n'
    random100.write(header)
    for file in os.listdir(Arc_dir):
        if file.split('.')[-1] == 'log':
            logfile = file
            #now we retrieve the best model chosen by AIC and BIC
            number_of_file =file.split('.fasta.log')[0].split('_')[-1]
            AIC_best_model = AIC_dic[number_of_file]
                    #now we retrieve the best model chosen by AIC and BIC
            BIC_best_model = BIC_dic[number_of_file]
            with open('{}/{}'.format(Arc_dir, logfile), 'r') as logout:
                for row in logout:
                    #check for the line with our best model
                    if AIC_best_model in row:
                        #this is to make sure that we are in model finder
                        if row[1].isnumeric():
                            list_row = row.split(' ')
                            #there are multiple blank spaces in this
                            for entry in list_row:
                                if entry == '':
                                    list_row.remove(entry)
                            for entry in list_row:
                                if entry == '':
                                    list_row.remove(entry)
                            
                            df_A = list_row[3]
                            AIC_score = list_row[4]
                            
                        #modelfinder has now started
                    #repeat same for BIC
                    if BIC_best_model in row:
                        #this is to make sure that we are in model finder
                        if row[1].isnumeric():
                            list_row = row.split(' ')
                            #there are multiple blank spaces in this
                            for entry in list_row:
                                if entry == '':
                                    list_row.remove(entry)
                            for entry in list_row:
                                if entry == '':
                                    list_row.remove(entry)
                            df_B = list_row[3]
                            BIC_score = list_row[-1].split('\n')[0]
            count = count + 1
            outline = number_of_file + '\t' +  AIC_best_model + '\t' +  AIC_score + '\t' +  df_A + '\t' +  BIC_best_model + '\t' +  BIC_score + '\t' +  df_B + '\n'
            Random100.write(outline)
                        #modelfinder has now started
