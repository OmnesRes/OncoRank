import scipy.stats as stats
import numpy as np
import math
from rpy2 import robjects as ro
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

cancers=['BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV',\
       'PAAD','READ','SARC','SKCM','STAD','UCEC']

genes={}

for cancer in cancers:
    ##path to https://github.com/OmnesRes/onco_lnc repository
    f=open(os.path.join(BASE_DIR,'onco_lnc','mrna','cox',cancer,'coeffs_pvalues.txt'))
    data=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
    for index,i in enumerate(data):
        genes[i[0]]=genes.get(i[0],[])+[(index+1)/float(len(data))-.5*1/float(len(data))]


chi_stats={}
for i in genes:
    stat=-2*sum([math.log(k) for k in genes[i]])
    chi_stats[i]=[len(genes[i]),stats.chi2.sf(stat,len(genes[i])*2)]

print len(chi_stats)



##script for creating a histogram

import pylab as plt

pvalues=[chi_stats[i][-1] for i in chi_stats if chi_stats[i][0]>=8]

##decide how man bins
number=100.0

counts={}

##use a dictionary to populate the bins
for i in range(int(number)):
    for j in pvalues:
        if i/number<j<=(i+1)/number:
            counts[i]=counts.get(i,0)+1

##convert the dictionary to a list
mylist=zip(counts.keys(),counts.values())

##sort the list so that the bins are in order
mylist.sort()


#plot the data with pylab
fig=plt.figure(figsize=(20, 8.3844))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=.21)
fig.subplots_adjust(left=.12)
fig.subplots_adjust(right=.98)


ax.bar([i[0]/number for i in mylist],[i[1] for i in mylist],color='b',width=1/number,linewidth=2.0)
ax.tick_params(axis='x',length=15,width=3,direction='out',labelsize=30)
ax.tick_params(axis='y',length=15,width=3,direction='out',labelsize=30)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['bottom'].set_position(['outward',10])
ax.spines['left'].set_position(['outward',10])
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xticks([i/10.0 for i in range(0,11)])
ax.set_xticklabels(['0']+[str(i/10.0) for i in range(1,11)])
ax.set_ylabel('Frequency',fontsize=50,labelpad=20)
ax.set_xlabel('Raw chi-squared p-value',fontsize=50,labelpad=20)
ax.set_xlim(-.005,1.005)
ax.spines['bottom'].set_bounds(0,1)
ax.set_title('mRNAs',fontsize=60,y=.9)
plt.savefig('mrnas.pdf')
plt.show()
