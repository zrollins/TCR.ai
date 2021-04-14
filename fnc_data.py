#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:24:58 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/fnc.csv', index_col = 0, header=None)
MART1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/fnc.csv', index_col = 0, header=None)
GVA = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/fnc.csv', index_col = 0, header=None)

#P2 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/2P/fnc.csv', index_col = 0, header=None)
#F3 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/3F/fnc.csv', index_col = 0, header=None)
#Y4 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/4Y/fnc.csv', index_col = 0, header=None)
#D5 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5D/fnc.csv', index_col = 0, header=None)
#H5 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5H/fnc.csv', index_col = 0, header=None)
#H6 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6H/fnc.csv', index_col = 0, header=None)
#Y6 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6Y/fnc.csv', index_col = 0, header=None)

Mtub2 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub2/fnc.csv', index_col = 0, header=None)
Mtub1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub1/fnc.csv', index_col = 0, header=None)
HSV1gp3 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/HSV1gp3/fnc.csv', index_col = 0, header=None)
S8 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/fnc.csv', index_col = 0, header=None)
V6 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/6V/fnc.csv', index_col = 0, header=None)
ImrA = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/ImrA/fnc.csv', index_col = 0, header=None)
hCD9 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/fnc.csv', index_col = 0, header=None)

data = data.append(MART1)
data = data.append(GVA)
#data = data.append(P2)
#data = data.append(F3)
#data = data.append(Y4)
#data = data.append(D5)
#data = data.append(H5)
#data = data.append(H6)
#data = data.append(Y6)
data = data.append(Mtub2)
data = data.append(Mtub1)
data = data.append(HSV1gp3)
data = data.append(S8)
data = data.append(V6)
data = data.append(ImrA)
data = data.append(hCD9)

print(data)

#Paramter weights
hb = data.iloc[0][141]
rxnd = data.iloc[0][142]
SASA = data.iloc[0][143]

w_rxnd = hb/rxnd
w_SASA = hb/SASA

data['Bond Strength'] = data[141]+ (w_rxnd*data[142]) + (w_SASA*data[143])
data = data.drop(labels=[141,142,143], axis=1)

print(data)
data.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/fnc_mstr_bio.csv', header=None)