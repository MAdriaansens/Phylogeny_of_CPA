#code for downloading all Euk proteomes expect for the ones from Joint Genome Institute or University of Ghent

#input Eukarya_metadata.tsv
#output below 245 proteomes downloaded, some need manual downloads, this is noted in the .tsv file itself. 

Euk_db = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_metadata.tsv'
import subprocess
JGI_list = []
with open(Euk_db, 'r') as E:
    for line in E:
        if line.split('\t')[-4] == 'JGI':
            JGI_list.append(line.split('\t')[1])
        else:
            print('direct downloading: {}'.format(line.split('\t')[1]))
            download_link = line.split('\t')[-1].split('\n')[0]
            subprocess.run(['wget', '{}'.format(download_link)])
print(JGI_list)
#Proteomonas_sulcata and Mastigamoeba_balamuthi need some small edits or help
