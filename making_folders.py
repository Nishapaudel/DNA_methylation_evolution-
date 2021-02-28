
        

import os
ly = [] 

with open('birds.tsv') as f:
    header = f.readline()
    for line in f:
        fields = line.strip( ).split('\t')
        new_name = fields[1].replace(' ','_')
        ly.append(new_name)

      

root_path = '/home/ccmb/Desktop/run_birds/birds_collection'
folders = ly
gh = folders


for folder in gh:
    try:
        os.mkdir(os.path.join(root_path,folder))
    except FileExistsError:
        pass
    
