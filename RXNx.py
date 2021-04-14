#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:09:49 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
##DO THIS FOR ALL STARTING CONFIGURATIONS##

t, t1, t2 = [], [], []
x, x1 = [], []
pep_c = []
pep_lj = []
MHC_c = []
MHC_lj = []

#COM Distance vs Time


with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            x1.append(float(cols[1]))
x = x1[::10]

#TCR-MHC IE vs Time

with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/ie_DEA_140ns.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            MHC_c.append(float(cols[1]))
            MHC_lj.append(float(cols[2]))
            
#TCR-peptide IE vs Time

with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/ie_DEC_140ns.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            pep_c.append(float(cols[1]))
            pep_lj.append(float(cols[2]))
            
ct = len(t) - len(x)            

print (len (x))
print(len (t1))
print(len (t))
print(len (t2))
print (len (MHC_c))
print (len (pep_c))

#PLOT Hbonds vs COM Distance

MHC_c = MHC_c[:len(MHC_c)-ct]
MHC_lj = MHC_lj[:len(MHC_lj)-ct]
pep_c = pep_c[:len(pep_c)-ct]
pep_lj = pep_lj[:len(pep_lj)-ct]

iex = {'COM Distance':x,'MHC-Coulombic':MHC_c,'MHC-LennardJones':MHC_lj,'Peptide_Coulombic':pep_c,'Peptide_LennardJones':pep_lj}

df=pd.DataFrame(iex)
df.to_excel('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/iex.xlsx')

