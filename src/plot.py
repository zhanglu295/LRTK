import matplotlib.pyplot as plt 
from collections import defaultdict
import os, sys, time

if len(sys.argv) != 3:
	sys.stderr.write("python %s fragment_info.csv outdir\n" % sys.argv[0])

csvfile = sys.argv[1]
pngdir = sys.argv[2]

#infile=open('fragment_info.csv')
infile=open(csvfile, 'r')
barcode_info=defaultdict(list)
Fragment_barcode=[]
Cr_list=[]
molecule_len=[]
for line in infile:
    A=line.strip('\n').split(',')
    barcode_info[A[1]].append(1)
    Cr_list.append(float(A[8]))
    molecule_len.append(int(A[5])/1000)
for k,v in barcode_info.items():
    if sum(v)>1:
        Fragment_barcode.append(sum(v))
#plot NFP
fig = plt.figure(figsize=(10,5))
plt.ylabel('Number of partitions',fontsize=20)
plt.xlabel('Fragments per partition',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.hist(barcode_info, 100, range=[0, 70], color = "skyblue")
y=[0, 20000,40000,60000,80000]
plt.yticks(y,['0','20,000','40,000','60,000','80,000']) 
plt.title('Distribution of the number of fragments per partition',fontsize=20)
plt.show()
NFP_png = os.path.join(pngdir, 'NFP.png')
#fig.savefig('NFP.png')
fig.savefig(NFP_png)

#plot Cr
import matplotlib.pyplot as plt 
fig = plt.figure(figsize=(10,5))
plt.ylabel('Number of fragments(M)',fontsize=20)
plt.xlabel('Depth',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.hist(Cr_list, 100, range=[0, 1.5], color = "skyblue")
y=[0, 200000,400000,600000,800000,1000000,1200000]
plt.yticks(y,['0','0.2','0.4','0.6','0.8','1','1.2']) 
plt.title('Distribution of sequencing depth per fragment',fontsize=20)
plt.show()
Cr_png = os.path.join(pngdir, 'Cr.png')
#fig.savefig('Cr.png')
fig.savefig(Cr_png)

#unweighted MuFL
fig = plt.figure(figsize=(10,5))
plt.ylabel('Number of fragments',fontsize=20)
plt.xlabel('Average fragment length(kbp)',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.hist(molecule_len, 100, range=[0, max(molecule_len)], color = "skyblue")
plt.title('Unweighted fragment length distribution',fontsize=20)
plt.show()
Unweight_MuFL_png = os.path.join(pngdir, 'Unweight_MuFL.png')
#fig.savefig('Unweight_MuFL.png')
fig.savefig(Unweight_MuFL_png)

#weighted MuFL
fig = plt.figure(figsize=(10,5))
plt.ylabel('Physical coverage(bp)',fontsize=20)
plt.xlabel('Average fragment length(kbp)',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.hist(molecule_len, 100, range=[0, max(molecule_len)],weights=molecule_len,color = "skyblue")
plt.title('Weighted fragment length distribution',fontsize=20)
plt.show()
Weight_MuFL_png = os.path.join(pngdir, 'Weight_MuFL.png')
#fig.savefig('Weight_MuFL.png')
fig.savefig(Weight_MuFL_png)
