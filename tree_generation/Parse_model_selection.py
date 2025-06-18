Euk_dir='/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/random_list'

Arc_dir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/intermediate_PF00999/random_list'

Bac_dir ='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/MERGED/random_list'

Seeds='/nesi/nobackup/uc04105/fasta_files/Seed_NhaA/randomList'

likelihood ='Seed_likelihood.txt'
model='Seeds_models.txt'

#each text file was generated using the grep 'Best-fit model:' subset_*/* and grep 'BEST SCORE FOUND' subset_*/* commands

model_dict = {}
score_dict = {} 

with open('{}/{}'.format(Seeds, model), 'r') as Mod:
    for line in Mod:
        numb =int((line.split('_')[1].split('/')[0]))
        model_chosen = (line.split(' chosen')[0].split(': ')[1])
        model_dict[numb] = model_chosen
with open('{}/{}'.format(Seeds, likelihood), 'r') as Like:
    for line in Like:
        numb =int((line.split('_')[1].split('/')[0]))
        score =(line.split('\n')[0].split(': ')[1])
        score_dict[numb] = score

total_dict = {}
for x in range(1,101):
    tlist = []
    model = model_dict[x]
    score = score_dict[x]
    name = 'Seed_{}'.format(x)
    tlist.append(model)
    tlist.append(score)
    total_dict[name] = tlist
score_Adict = {}
model_Adict = {}

likelihoodA = 'Archaea_likelihood.txt'
modelA='Archaea_tree_model.txt'

with open('{}/{}'.format(Arc_dir, modelA), 'r') as Mod:
    for line in Mod:
        numb =int((line.split('_')[1].split('/')[0]))
        model_chosen = (line.split(' chosen')[0].split(': ')[1])
        model_Adict[numb] = model_chosen
with open('{}/{}'.format(Arc_dir, likelihoodA), 'r') as Like:
    for line in Like:
        numb =int((line.split('_')[1].split('/')[0]))
        score =(line.split('\n')[0].split(': ')[1])
        score_Adict[numb] = score

for x in range(1,101):
    tlist = []
    model = model_Adict[x]
    score = score_Adict[x]
    name = 'Arc_{}'.format(x)
    tlist.append(model)
    tlist.append(score)
    total_dict[name] = tlist
  score_Bdict = {}
model_Bdict = {}

likelihoodB ='Best_score.txt'
modelB= 'Bacteria_hmmscanned_model.txt'

with open('{}/{}'.format(Bac_dir, modelB), 'r') as Mod:
    for line in Mod:
        numb =int((line.split('_')[1].split('/')[0]))
        model_chosen = (line.split(' chosen')[0].split(': ')[1])
        model_Bdict[numb] = model_chosen
with open('{}/{}'.format(Bac_dir, likelihoodB), 'r') as Like:
    for line in Like:
        numb =int((line.split('_')[1].split('/')[0]))
        score =(line.split('\n')[0].split(': ')[1])
        score_Bdict[numb] = score

for x in range(1,101):
    tlist = []
    model = model_Bdict[x]
    score = score_Bdict[x]
    name = 'Bac_{}'.format(x)
    tlist.append(model)
    tlist.append(score)
    total_dict[name] = tlist
score_Edict = {}
model_Edict = {}

likelihoodE ='Eukarya_likelihood.txt'
modelE= 'Model_picked_euk_pf00999.txt'

with open('{}/{}'.format(Euk_dir, modelE), 'r') as Mod:
    for line in Mod:
        numb =int((line.split('_')[1].split('/')[0]))
        model_chosen = (line.split(' chosen')[0].split(': ')[1])
        model_Edict[numb] = model_chosen
with open('{}/{}'.format(Euk_dir, likelihoodE), 'r') as Like:
    for line in Like:
        numb =int((line.split('_')[1].split('/')[0]))
        score =(line.split('\n')[0].split(': ')[1])
        score_Edict[numb] = score

for x in range(1,101):
    tlist = []
    model = model_Edict[x]
    score = score_Edict[x]
    name = 'Euk_{}'.format(x)
    tlist.append(model)
    tlist.append(score)
    total_dict[name] = tlist


print(len(list(total_dict.keys())))
#should total 400. 
