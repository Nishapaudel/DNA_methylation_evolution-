
import random
import py2bit
import sys
import os


print('Chrom', 'Feature', 'start', 'end', 'CpG_ratio', 'GpC_ratio', sep = '\t')
bit_file =sys.argv[1]
tb = py2bit.open('bit_file', True)
tb_dict =tb.chroms()
chrom = list(tb_dict.keys()) #chrom_names
len_chrom = list(tb_dict.values()) #chrom_lens
count_entry = 0

while count_entry < 10000:
    j = random.choice(chrom)
    random_start = random.randrange(1,tb_dict[j])
    random_end = random_start + 1000
    if random_end > int(tb_dict[j]):
              continue
    else:
        perseq = tb.sequence(j, random_start, random_end)
        C = int(perseq.upper().count('C'))
        G = int(perseq.upper().count('G'))
        CG = int(perseq.upper().count('CG'))
        GC = int(perseq.upper().count('GC'))
        if GC < 1 or CG < 1 or C < 1 or G < 1:
            continue
        CpG_ratio = round((CG/1000)/((C/1000)*(G/1000)),4)
        GpC_ratio = round((GC/1000)/((C/1000)*(G/1000)),4)
        print(j,'random seq' ,random_start, random_end, CpG_ratio, GpC_ratio, sep = '\t')
        count_entry = count_entry + 1
        
