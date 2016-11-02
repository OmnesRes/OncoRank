##script for finding the overlap in the top 5% most significant genes in each cancer and plotting results

##load necessary modules
import pylab as plt
import numpy as np
import math
import os
from itertools import *


##BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
##
####get the 100 most significant genes for each cancer
##
####path to https://github.com/OmnesRes/onco_lnc repository
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','BLCA','coeffs_pvalues.txt'))
##BLCA=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','BRCA','coeffs_pvalues.txt'))
##BRCA=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##                           
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','CESC','coeffs_pvalues.txt'))
##CESC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','COAD','coeffs_pvalues.txt'))
##COAD=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','ESCA','coeffs_pvalues.txt'))
##ESCA=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','HNSC','coeffs_pvalues.txt'))
##HNSC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','KIRC','coeffs_pvalues.txt'))
##KIRC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##                           
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','KIRP','coeffs_pvalues.txt'))
##KIRP=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','LAML','coeffs_pvalues.txt'))
##LAML=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','LGG','coeffs_pvalues.txt'))
##LGG=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','LIHC','coeffs_pvalues.txt'))
##LIHC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','LUAD','coeffs_pvalues.txt'))
##LUAD=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','LUSC','coeffs_pvalues.txt'))
##LUSC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','OV','coeffs_pvalues.txt'))
##OV=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','PAAD','coeffs_pvalues.txt'))
##PAAD=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','READ','coeffs_pvalues.txt'))
##READ=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','SARC','coeffs_pvalues.txt'))
##SARC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','SKCM','coeffs_pvalues.txt'))
##SKCM=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','STAD','coeffs_pvalues.txt'))
##STAD=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
##f=open(os.path.join(BASE_DIR,'onco_lnc','mirna','cox','UCEC','coeffs_pvalues.txt'))
##UCEC=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
##
###determine genes in the population and make dictionaries
##
##all_cancers=[BLCA,BRCA,CESC,COAD,ESCA,HNSC,KIRC,KIRP,LAML,LGG,LIHC,LUAD,LUSC,OV,\
##             PAAD,READ,SARC,SKCM,STAD,UCEC]
##
##
##population={}
##for i in all_cancers:
##    for j in i:
##        population[j[0]]=''
##
##print len(population)
####figure out 5th percentile
##per_5=int(.05*len(population))
##print per_5
##values=[]
##
##for size in range(2,6):
##    print size
##    temp=[]
##    count=0
##    for i in combinations(all_cancers,size):
##        combined_list=[]
##        for j in i:
##            combined_list+=j[:per_5]
##        counts={}
##        for k in combined_list:
##            counts[k[0]]=counts.get(k[0],0)+1
##        count+=1
##        temp.append(len([gene for gene in counts if counts[gene]==size]))
##    values.append(temp)
##
##f=open('mirna_values.txt','w')
##f.write(str(values))
##f.close()


#5_per=34

f=open('mirna_values.txt')
values=eval(f.read())
values=[list(np.array(i)/34.0*100) for i in values]
fig=plt.figure(figsize=(20, 8.3844))  
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=.18)
fig.subplots_adjust(left=.07)
fig.subplots_adjust(right=1.00)

positions=[(i)/2.0 for i in range(len(values))]
violin_parts=ax.violinplot(values,positions=positions,showextrema=False,widths=.4)
for pc in violin_parts['bodies']:
    pc.set_facecolor('red')
    pc.set_edgecolor('black')
    pc.set_linewidths(1.5)
    pc.set_alpha(1.0)



ax.set_xticks(positions)
ax.set_xticklabels([str(i+2) for i in range(len(values))])
ax.tick_params(axis='y',length=0,width=0,direction='out',labelsize=30)
ax.tick_params(axis='x',length=0,labelsize=40,pad=15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_position(['outward',0])
ax.yaxis.set_ticks_position('left')
ax.set_ylabel('Genes shared (%)',fontsize=50,labelpad=20)
ax.set_xlabel('Cancers compared',fontsize=50,labelpad=20)
ax.set_xlim(-.25,1.75)
ax.set_ylim(-1,50)
ax.set_title('miRNAs',fontsize=60,y=.9)
plt.savefig('mirnas.pdf')
plt.show()


