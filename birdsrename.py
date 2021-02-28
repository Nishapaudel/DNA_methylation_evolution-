import os
names = {}
with open('birds.tsv') as f:
    header= f.readline()
    for line in f:
        fields = line.strip().split('\t')
        new_name = fields[1].replace(' ','_')+'.gff'
        old_name = fields[-1] + '.gff'
        names[old_name] = new_name
for n in names:
   try:
       os.rename(n, names[n])
   except FileNotFoundError:
       pass
