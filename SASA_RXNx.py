#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 10:01:16 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#Residues

resi = ['Y','S','D','R','G','S','Q','S','F',
        'Y','S','N','G','D','K',
        'A','V','N','F','G','G','G','K','L','I','F',
        'Q','D','M','R','H','N','A',
        'Y','S','N','T','A','G','T','T',
        'A','S','S','L','S','F','G','T','E','A','F','F',
        'M','A','A','Q','T','T','K','H','K','W','E','A',
        'A','H','V','A','E','Q','L','R','A','Y','L','E',
        'G','T','C','V','E','W','L','R','R','Y','L','E',
        'N','G','K','E','T','L',
        'P','W','I','E','Q','E','G','P','E','Y','W','D',
        'G','E','T','R','K','V','K','A','H','S','Q','T',
        'H','R','V','D','L','G','T','L','R','G','Y','Y',
        'A','V','G','I','G','I','A','V','V']

#SASA
t, t1, t2, pep, pep1, pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/SASA_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/130/SASA_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/120/SASA_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            pep2.append(float(cols[1]))

#Bin SASA  
cut_bins = [12.0,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13.0,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14.0,14.1,14.2]
            
sasa={'pep (1ns)': pep}
df = pd.DataFrame(sasa)
s=pd.cut(df['pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'pep (2ns)': pep1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'pep (3ns)': pep2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['pep (2ns)'] = s1['pep (2ns)']
s2 = s2.to_frame()
SASA['pep (3ns)'] = s2['pep (3ns)']

l = [12.0,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13.0,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14.0,14.1]
x = pd.Series(l)
SASA['SASA']= x.values

M1 = ['pep (1ns)','pep (2ns)','pep (3ns)']
SASA['pep'] = SASA[M1].astype(float).mean(axis=1)
SASA['pep SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx = SASA['pep'].idxmax().right

print(mx)

#RXN Distance

df = pd.read_excel (r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/iex.xlsx')
df2 = pd.read_excel (r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/130/iex.xlsx')
df3 = pd.read_excel (r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/120/iex.xlsx')

df = df.append(df2, ignore_index=True)
df = df.append(df3, ignore_index=True)

bins=np.linspace(5.884,14.612,200).tolist()
MHCc = df.groupby(pd.cut(df['COM Distance'], bins=bins))['MHC-Coulombic'].agg(['mean','sem','size'])
MHClj = df.groupby(pd.cut(df['COM Distance'], bins=bins))['MHC-LennardJones'].agg(['mean','sem','size'])
pepc = df.groupby(pd.cut(df['COM Distance'], bins=bins))['Peptide_Coulombic'].agg(['mean','sem','size'])
peplj = df.groupby(pd.cut(df['COM Distance'], bins=bins))['Peptide_LennardJones'].agg(['mean','sem','size'])


zero_MHCc = MHCc['mean'].idxmax().right
zero_MHClj = MHClj['mean'].idxmax().right
zero_pepc = pepc['mean'].idxmax().right
zero_peplj = peplj['mean'].idxmax().right

zero = max(zero_MHCc,zero_MHClj,zero_pepc,zero_peplj)

# Total H-Bonds

hmap1 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/hmap.csv', index_col = 0, header=None)
hmap2 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/130/hmap.csv', index_col = 0, header=None)
hmap3 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/120/hmap.csv', index_col = 0, header=None)

hb1, t1 = hmap1.shape
hb2, t2 = hmap2.shape
hb3, t3 =  hmap3.shape

hb = (hb1 + hb2 + hb3)/3

#Compile Data for RNN

datar_labels = {'hCD9':resi}
datar = pd.DataFrame(datar_labels)
print(datar)
datar = datar.replace(['A','R','N','D','C','Q','E','G','H','I','L',
                         'K','M','F','P','S','T','W','Y','V']
                       ,['1','2','3','4','5','6','7','8','9','10',
                       '11','12','13','14','15','16','17','18','19','20'])

datar = datar.T
datar['Total H-Bonds'] = hb
datar['RXN Distance'] = zero
datar['SASA'] = mx


datar.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/fnc.csv', header=None)

