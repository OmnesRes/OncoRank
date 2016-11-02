import scipy.stats as stats
import numpy as np
import math
from rpy2 import robjects as ro
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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


genes_reversed={}
for cancer in cancers:
    ##path to https://github.com/OmnesRes/onco_lnc repository
    f=open(os.path.join(BASE_DIR,'onco_lnc','mrna','cox',cancer,'coeffs_pvalues.txt'))
    data=sorted([[j.split('\t')[0],float(j.split('\t')[-2])] for j in f],key=lambda x:x[-1])
    for index,i in enumerate(data):
        genes_reversed[i[0]]=genes_reversed.get(i[0],[])+[(len(data)-index)/float(len(data))-.5*1/float(len(data))]


chi_stats_reversed={}
for i in genes_reversed:
    stat=-2*sum([math.log(k) for k in genes_reversed[i]])
    chi_stats_reversed[i]=[len(genes_reversed[i]),stats.chi2.sf(stat,len(genes_reversed[i])*2)]



print len(chi_stats)
print len(chi_stats_reversed)

merged=[]
for i in chi_stats:
    if chi_stats[i][0]>=8:
        if chi_stats[i][1]<chi_stats_reversed[i][1]:
            merged.append([i,math.log(2*chi_stats[i][1],10),2*chi_stats[i][1]])
        elif chi_stats[i][1]>chi_stats_reversed[i][1]:
            merged.append([i,-1*math.log(2*chi_stats_reversed[i][1],10),2*chi_stats_reversed[i][1]])
        else:
            merged.append([i,0,2*chi_stats_reversed[i][1]])
            print i,chi_stats_reversed[i][1]
            
print len(merged)


f=open(os.path.join(BASE_DIR,'onco_rank','oncolnc_mrna_names.txt'))
data=[eval(i.strip()) for i in f]
name_dict={}
for i in data:
    name_dict[i[0]]=i[1:]

print len(name_dict)

merged.sort(key=lambda x:x[-1])
##BH correction
pvalues=[i[-1] for i in merged]
#prepare data for R
ro.globalenv['pvalues']=ro.FloatVector(pvalues)
#perform BH adjustment
res=ro.r('p.adjust(pvalues,"BH")')
#extract data
adjusted=list(res)


f=open('mrna_ranks.txt','w')
for i,fdr in zip(merged,adjusted):
    f.write(i[0])
    f.write('\t')
    for j in name_dict[i[0]]:
        f.write(j)
        f.write('\t')
    f.write(str(i[1]))
    f.write('\t')
    f.write(str(fdr))
    f.write('\n')
f.close()
