
import os
import sys 
import py2bit
import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import dim
hv.extension('bokeh')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
#====================================File 1: generates the ratio==================================================#
def extract_promos(gff_file, bd):

    ''' 
    Get the coordinates from gff file, fetch the sequence from 2bit file and calculate the ratio 
    '''
    #output file
    fout_promo = open ('promoters_output.txt' , 'w')
    #print headings
    print('Chromosome', 'Feature', 'Start', 'End', 'Strand', 'CpG_ratio', 'GpC_ratio', sep = '\t', file = fout_promo)

    with open(gff_file) as fh:

        chrom_name = ''
        chrom_seq_length = 0
   
        for line in fh:
            if line.startswith('#'):
                    continue
            content = line.strip().split('\t')
            if chrom_name == '' or chrom_name != content[0]:
                chrom_name = content[0]
                chrom_seq_length = len(bd.sequence(str(content[0])))
          
            chrom= content[0]
            feature = content[2]
            start = int(content[3])
            end = int(content[4])
            strand = content[6]
            window_end = 0
            window_start = 0

            if feature == 'gene':
                if strand == '+':
                    window_start = start - 1000
                    window_end = window_start + 2000
                    if window_start > window_end or window_start < 1 or window_end < 1 or window_start == window_end or window_start > chrom_seq_length or window_end > chrom_seq_length or (window_end - window_start) < 1:
                        continue
                    else:
                        perseq = bd.sequence(chrom,window_start, window_end)
                        

                else:                
                    window_start = end - 1000
                    window_end = window_start + 2000
                    if window_start > window_end or window_start < 1 or window_end < 1 or window_start == window_end or window_start > chrom_seq_length or window_end > chrom_seq_length or (window_end-window_start) < 1:
                        continue
                        
                    else:
                        perseq = bd.sequence(chrom,window_start, window_end)
                        C = int(perseq.upper().count('C'))
                        G = int(perseq.upper().count('G'))
                        CG = int(perseq.upper().count('CG'))
                        GC = int(perseq.upper().count('GC'))
                        seq_len = len(perseq)
                        if GC < 1 or CG < 1 or G < 1 or C < 1 or seq_len < 1:
                            continue
                        
                        CpG_ratio = round(((CG / (G*C))*seq_len),2)
                        GpC_ratio = round(((GC / (G*C))*seq_len),2)
                print(chrom, feature, start, end,  strand, CpG_ratio, GpC_ratio, sep = '\t', file = fout_promo)

#main body
bit_file = sys.argv[1]
gff_file = sys.argv[2]

bd = py2bit.open(bit_file , True)
extract_promos(gff_file, bd)

#closing the files
bd.close()

#=============================file 2 = arranged data for violin plot============================================#

out_file = 'promoters_output.txt'

# arraning work done ..

def violin_plot(out_file):
    Violin_out = open('arranged_data_violin.txt', 'w')
    print('Dinucleotide', 'ratio', sep= '\t', file = Violin_out)
    with open (out_file) as fh:
        header = fh.readline()
        
        for line in fh:
            columns = line.strip().split('\t')
            CpG =float(columns[5])
            print( 'CpG' , CpG, sep = '\t', file = Violin_out)
            GpC =float(columns[6])
            print( 'GpC' , GpC, sep = '\t', file = Violin_out)
            

violin_plot(out_file)

#==============================plot 1 = violin plot==============================


violin_data = 'arranged_data_violin.txt'

#loading the data
df = pd.read_csv(violin_data, sep = '\t')
# Declearing the data
#plot
violin = hv.Violin(df, ('Dinucleotide', 'Dinucleotides'), ('ratio', 'Ratio')).redim.range(ratio=(0, 2))
violin.opts(height=500, width=900, violin_fill_color=dim('Dinucleotides').str(), cmap='Set1')
#save plot
hv.save(violin, 'Promo_violin_plot.png', fmt = 'png' )

#========================making the line plot with guassian filter=========================


count_CpG = {}
count_GpC = {}

with open ('promoters_output.txt') as fh:
    header = fh.readline()
    for line in fh:
        content = line.strip().split('\t')
        CpG_ratio = float(content[5])
        GpC_ratio = float(content[6])
        if CpG_ratio > 2:
            continue
       
        if CpG_ratio in count_CpG.keys():
            count_CpG[CpG_ratio] += 1
        else:
            count_CpG[CpG_ratio] = 1
       
        
        if GpC_ratio > 2:
            continue
       
        if GpC_ratio in count_GpC.keys():
            count_GpC[GpC_ratio] += 1
        else:
            count_GpC[GpC_ratio] = 1
       
plt.figure(figsize = (10,7))
plt.title('Methylation')
x1 = list(count_CpG.keys())
x1.sort()
x1.pop(0)
#print(x1)
y1 = []
for i in x1:
    y1.append(count_CpG[i])
x2 = list(count_GpC.keys())
x2.sort()
x2.pop(0)
y2 = []
for i in x2:
    y2.append(count_GpC[i])
#y1 = [i for i in count_CpG.values()]
#print(list(y1))
guss_y1 = gaussian_filter1d(y1, sigma = 3)
plt.fill_between(x1, y1, alpha = 0.3, color = 'c')
plt.plot(x1, guss_y1, color = 'c', label = 'CpG')
peak = find_peaks(guss_y1)

guss_y2 = gaussian_filter1d(y2, sigma = 3)
plt.fill_between(x2, y2, alpha = 0.3, color = 'y')
plt.plot(x2, guss_y2, color = 'y', label = 'GpC')
peak = find_peaks(guss_y2)

for i in peak[0]:
    plt.plot(x1[i], guss_y1[i], '|', color = 'black', markersize = 15)
for i in peak[0]:
    plt.plot(x2[i], guss_y2[i], '|', color = 'black', markersize = 15)
plt.xlabel('ratio')
plt.ylabel('frequency')
plt.legend()
plt.savefig('promo_lineplot.png')

#====================================making histograms========================



#====================================skewness and kurtosis====================

with open('skew_kurt.txt','w') as fout_skew:
    #Read the data set
    data=pd.read_csv("promoters_output.txt", sep = "\t")
    #creating data frame
    df=pd.DataFrame(data)

    print( "The value of Skewness is:", file =fout_skew)
    #calculating the skewness
    print(df.skew(), file = fout_skew)
    print("The value of kurtosis is:", file = fout_skew)
    #calculating the kurtosis
    print(df.kurtosis(), file= fout_skew)



