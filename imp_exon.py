
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



def extract_exons(gff_file, bd):
    ''' 
    Get the coordinates from gff file, fetch the sequence from 2bit file and calculate the ratio 
    '''
    
    fout_exons = open('exons_output.txt', 'w')
    print("Chromosome","Feature","Start","End","Strand", "CpG_ratio", "Gpc_ratio" , sep='\t', file = fout_exons)
    with open(gff_file) as fh:
        chr_name = ''
        chr_seq_length = 0
        for line in fh:
            if line.startswith('#'):
                continue
            content = line.strip().split('\t')
            if chr_name == '' or chr_name != content[0]:
                chr_name = content[0]
                chr_seq_length = len(bd.sequence(str(content[0])))
            chrom = content[0]
            feature = content[2]
            start = int(content[3])
            end = int(content[4])
       
       
            
            if feature == 'exon':
                exon_start = start
                exon_end = end
                if exon_start > exon_end or exon_start < 1 or exon_end < 1 or exon_start == exon_end or (exon_end-exon_start) < 1 or exon_start > chr_seq_length or exon_end > chr_seq_length :
                    continue
                perseq = bd.sequence(chrom,exon_start,exon_end)
                
                G = int(perseq.upper().count('G'))
                C = int(perseq.upper().count('C'))
               
          
                GC = int(perseq.upper().count('GC'))
               
                CG = int(perseq.upper().count('CG'))
              
          
                seq_len = len(perseq)
                if GC < 1 or  CG < 1 or   G < 1 or C < 1 or seq_len < 100:
                    continue
               
                GC_ratio = round(((GC/ (G*C))*seq_len),2)
 
                CG_ratio = round(((CG/ (C*G))*seq_len),2)
               
                print(content[0],feature, exon_start, exon_end,content[6], CG_ratio, GC_ratio , sep = '\t', file = fout_exons )
        fout_exons.close()
              

#main body
bit_file = sys.argv[1]
gff_file = sys.argv[2]
bd = py2bit.open(bit_file, True)
extract_exons(gff_file, bd)
bd.close()



#=============================file 2 = arranged data for violin plot============================================#

out_file = 'exons_output.txt'

# arraning work done ..

def violin_plot(out_file):
    Violin_out = open('exon_arranged_data_violin.txt', 'w')
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


violin_data = 'exon_arranged_data_violin.txt'

#loading the data
df = pd.read_csv(violin_data, sep = '\t')
# Declearing the data
#plot
violin = hv.Violin(df, ('Dinucleotide', 'Dinucleotides'), ('ratio', 'Ratio')).redim.range(ratio=(0, 2))
violin.opts(height=500, width=900, violin_fill_color=dim('Dinucleotides').str(), cmap='Set1')
#save plot
hv.save(violin, 'exon_violin_plot.png', fmt = 'png' )

#========================making the line plot with guassian filter=========================


count_CpG = {}
count_GpC = {}

with open ('exons_output.txt') as fh:
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
plt.savefig('exon_lineplot.png')

#====================================making histograms========================



#====================================skewness and kurtosis====================

with open('exon_skew_kurt.txt','w') as fout_skew:
    #Read the data set
    data=pd.read_csv("expns_output.txt", sep = "\t")
    #creating data frame
    df=pd.DataFrame(data)

    print( "The value of Skewness is:", file =fout_skew)
    #calculating the skewness
    print(df.skew(), file = fout_skew)
    print("The value of kurtosis is:", file = fout_skew)
    #calculating the kurtosis
    print(df.kurtosis(), file= fout_skew)



