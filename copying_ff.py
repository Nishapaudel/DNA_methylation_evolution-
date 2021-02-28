import shutil,glob
files = glob.glob('/home/ccmb/Desktop/run_birds/birdss_gff/*.gff')
files.sort()
folders = glob.glob('/home/ccmb/Desktop/run_birds/birds_collection/*/')
folders.sort()
for i in range(len(files)):
    if files[i].split('/')[-1].split('.')[0] == folders[i].split('/') [-2]:
        shutil.copy(files[i],folders[i])
