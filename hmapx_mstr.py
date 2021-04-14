#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 21:35:56 2021

@author: zrollins
"""
import numpy as np
import pandas as pd
import statistics as st
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import gromacs
from gromacs.formats import XPM

L1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/hmapx_pmhc.csv', index_col = 0, header=None)
MART1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/hmapx_pmhc.csv', index_col = 0, header=None)
GVA = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/hmapx_pmhc.csv', index_col = 0, header=None)

P2 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/2P/hmapx_pmhc.csv', index_col = 0, header=None)
F3 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/3F/hmapx_pmhc.csv', index_col = 0, header=None)
Y4 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/4Y/hmapx_pmhc.csv', index_col = 0, header=None)
D5 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5D/hmapx_pmhc.csv', index_col = 0, header=None)
H5 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5H/hmapx_pmhc.csv', index_col = 0, header=None)
H6 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6H/hmapx_pmhc.csv', index_col = 0, header=None)
Y6 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6Y/hmapx_pmhc.csv', index_col = 0, header=None)

Mtub2 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub2/hmapx_pmhc.csv', index_col = 0, header=None)
ImrA = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/ImrA/hmapx_pmhc.csv', index_col = 0, header=None)
hCD9 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/hmapx_pmhc.csv', index_col = 0, header=None)
Mtub1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub1/hmapx_pmhc.csv', index_col = 0, header=None)
HSV1gp3 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/HSV1gp3/hmapx_pmhc.csv', index_col = 0, header=None)
S8 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/hmapx_pmhc.csv', index_col = 0, header=None)
V6 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/6V/hmapx_pmhc.csv', index_col = 0, header=None)

#Indices
ndx_L1 = pd.Index(L1.index.tolist())
ndx_M1 = pd.Index(MART1.index.tolist())
ndx_GVA = pd.Index(GVA.index.tolist())

ndx_P2 = pd.Index(P2.index.tolist())
ndx_F3 = pd.Index(F3.index.tolist())
ndx_Y4 = pd.Index(Y4.index.tolist())
ndx_D5 = pd.Index(D5.index.tolist())
ndx_H5 = pd.Index(H5.index.tolist())
ndx_H6 = pd.Index(H6.index.tolist())
ndx_Y6 = pd.Index(Y6.index.tolist())

ndx_Mtub2 = pd.Index(Mtub2.index.tolist())
ndx_ImrA = pd.Index(ImrA.index.tolist())
ndx_hCD9 = pd.Index(hCD9.index.tolist())
ndx_Mtub1 = pd.Index(Mtub1.index.tolist())
ndx_HSV1gp3 = pd.Index(HSV1gp3.index.tolist())
ndx_S8 = pd.Index(S8.index.tolist())
ndx_V6 = pd.Index(V6.index.tolist())



#Total Common Index
L1_M1  = ndx_M1.difference(ndx_L1, sort=False).tolist()
L1_GVA = ndx_GVA.difference(ndx_L1, sort=False).tolist()
L1_P2 = ndx_P2.difference(ndx_L1, sort=False).tolist()
L1_F3 = ndx_F3.difference(ndx_L1, sort=False).tolist()
L1_Y4 = ndx_Y4.difference(ndx_L1, sort=False).tolist()
L1_D5 = ndx_D5.difference(ndx_L1, sort=False).tolist()
L1_H5 = ndx_H5.difference(ndx_L1, sort=False).tolist()
L1_H6 = ndx_H6.difference(ndx_L1, sort=False).tolist()
L1_Y6 = ndx_Y6.difference(ndx_L1, sort=False).tolist()
L1_Mtub2 = ndx_Mtub2.difference(ndx_L1, sort=False).tolist()
L1_ImrA = ndx_ImrA.difference(ndx_L1, sort=False).tolist()
L1_hCD9 = ndx_hCD9.difference(ndx_L1, sort=False).tolist()
L1_Mtub1 = ndx_Mtub1.difference(ndx_L1, sort=False).tolist()
L1_HSV1gp3 = ndx_HSV1gp3.difference(ndx_L1, sort=False).tolist()
L1_S8 = ndx_S8.difference(ndx_L1, sort=False).tolist()
L1_V6 = ndx_V6.difference(ndx_L1, sort=False).tolist()
L1_diff= list(set(L1_M1 + L1_GVA + L1_P2 + L1_F3 + L1_Y4 + L1_D5 + L1_H5 + L1_H6 + L1_Y6 
                  + L1_Mtub2 + L1_Mtub1 + L1_HSV1gp3 + L1_S8 + L1_V6 + L1_hCD9 + L1_ImrA))

M1_L1 = M1_L1= ndx_L1.difference(ndx_M1, sort=False).tolist()
M1_GVA= ndx_GVA.difference(ndx_M1, sort=False).tolist()
M1_P2 = ndx_P2.difference(ndx_M1, sort=False).tolist()
M1_F3 = ndx_F3.difference(ndx_M1, sort=False).tolist()
M1_Y4 = ndx_Y4.difference(ndx_M1, sort=False).tolist()
M1_D5 = ndx_D5.difference(ndx_M1, sort=False).tolist()
M1_H5 = ndx_H5.difference(ndx_M1, sort=False).tolist()
M1_H6 = ndx_H6.difference(ndx_M1, sort=False).tolist()
M1_Y6 = ndx_Y6.difference(ndx_M1, sort=False).tolist()
M1_Mtub2 = ndx_Mtub2.difference(ndx_M1, sort=False).tolist()
M1_ImrA = ndx_ImrA.difference(ndx_M1, sort=False).tolist()
M1_hCD9 = ndx_hCD9.difference(ndx_M1, sort=False).tolist()
M1_Mtub1 = ndx_Mtub1.difference(ndx_M1, sort=False).tolist()
M1_HSV1gp3 = ndx_HSV1gp3.difference(ndx_M1, sort=False).tolist()
M1_S8 = ndx_S8.difference(ndx_M1, sort=False).tolist()
M1_V6 = ndx_V6.difference(ndx_M1, sort=False).tolist()
M1_diff= list(set(M1_L1 + M1_GVA + M1_P2 + M1_F3 + M1_Y4 + M1_D5 + M1_H5 + M1_H6 + M1_Y6 
                  + M1_Mtub2 + M1_Mtub1 + M1_HSV1gp3 + M1_S8 + M1_V6 + M1_ImrA + M1_hCD9))

GVA_M1= ndx_M1.difference(ndx_GVA, sort=False).tolist()
GVA_L1= ndx_L1.difference(ndx_GVA, sort=False).tolist()
GVA_P2 = ndx_P2.difference(ndx_GVA, sort=False).tolist()
GVA_F3 = ndx_F3.difference(ndx_GVA, sort=False).tolist()
GVA_Y4 = ndx_Y4.difference(ndx_GVA, sort=False).tolist()
GVA_D5 = ndx_D5.difference(ndx_GVA, sort=False).tolist()
GVA_H5 = ndx_H5.difference(ndx_GVA, sort=False).tolist()
GVA_H6 = ndx_H6.difference(ndx_GVA, sort=False).tolist()
GVA_Y6 = ndx_Y6.difference(ndx_GVA, sort=False).tolist()
GVA_Mtub2 = ndx_Mtub2.difference(ndx_GVA, sort=False).tolist()
GVA_ImrA = ndx_ImrA.difference(ndx_GVA, sort=False).tolist()
GVA_hCD9 = ndx_hCD9.difference(ndx_GVA, sort=False).tolist()
GVA_Mtub1 = ndx_Mtub1.difference(ndx_GVA, sort=False).tolist()
GVA_HSV1gp3 = ndx_HSV1gp3.difference(ndx_GVA, sort=False).tolist()
GVA_S8 = ndx_S8.difference(ndx_GVA, sort=False).tolist()
GVA_V6 = ndx_V6.difference(ndx_GVA, sort=False).tolist()
GVA_diff= list(set(GVA_M1 + GVA_L1 + GVA_P2 + GVA_F3 + GVA_Y4 + GVA_D5 + GVA_H5 + GVA_H6 + GVA_Y6
                   + GVA_Mtub2 + GVA_Mtub1 +  GVA_HSV1gp3 + GVA_S8 + GVA_V6 + GVA_ImrA + GVA_hCD9))

P2_M1= ndx_M1.difference(ndx_P2, sort=False).tolist()
P2_L1= ndx_L1.difference(ndx_P2, sort=False).tolist()
P2_GVA = ndx_GVA.difference(ndx_P2, sort=False).tolist()
P2_F3 = ndx_F3.difference(ndx_P2, sort=False).tolist()
P2_Y4 = ndx_Y4.difference(ndx_P2, sort=False).tolist()
P2_D5 = ndx_D5.difference(ndx_P2, sort=False).tolist()
P2_H5 = ndx_H5.difference(ndx_P2, sort=False).tolist()
P2_H6 = ndx_H6.difference(ndx_P2, sort=False).tolist()
P2_Y6 = ndx_Y6.difference(ndx_P2, sort=False).tolist()
P2_Mtub2 = ndx_Mtub2.difference(ndx_P2, sort=False).tolist()
P2_ImrA = ndx_ImrA.difference(ndx_P2, sort=False).tolist()
P2_hCD9 = ndx_hCD9.difference(ndx_P2, sort=False).tolist()
P2_Mtub1 = ndx_Mtub1.difference(ndx_P2, sort=False).tolist()
P2_HSV1gp3 = ndx_HSV1gp3.difference(ndx_P2, sort=False).tolist()
P2_S8 = ndx_S8.difference(ndx_P2, sort=False).tolist()
P2_V6 = ndx_V6.difference(ndx_P2, sort=False).tolist()
P2_diff= list(set(P2_M1 + P2_L1 + P2_GVA + P2_F3 + P2_Y4 + P2_D5 + P2_H5 + P2_H6 + P2_Y6
                  + P2_Mtub2 + P2_Mtub1 +  P2_HSV1gp3 + P2_S8 + P2_V6 + P2_ImrA + P2_hCD9))

F3_M1= ndx_M1.difference(ndx_F3, sort=False).tolist()
F3_L1= ndx_L1.difference(ndx_F3, sort=False).tolist()
F3_GVA = ndx_GVA.difference(ndx_F3, sort=False).tolist()
F3_P2 = ndx_P2.difference(ndx_F3, sort=False).tolist()
F3_Y4 = ndx_Y4.difference(ndx_F3, sort=False).tolist()
F3_D5 = ndx_D5.difference(ndx_F3, sort=False).tolist()
F3_H5 = ndx_H5.difference(ndx_F3, sort=False).tolist()
F3_H6 = ndx_H6.difference(ndx_F3, sort=False).tolist()
F3_Y6 = ndx_Y6.difference(ndx_F3, sort=False).tolist()
F3_Mtub2 = ndx_Mtub2.difference(ndx_F3, sort=False).tolist()
F3_ImrA = ndx_ImrA.difference(ndx_F3, sort=False).tolist()
F3_hCD9 = ndx_hCD9.difference(ndx_F3, sort=False).tolist()
F3_Mtub1 = ndx_Mtub1.difference(ndx_F3, sort=False).tolist()
F3_HSV1gp3 = ndx_HSV1gp3.difference(ndx_F3, sort=False).tolist()
F3_S8 = ndx_S8.difference(ndx_F3, sort=False).tolist()
F3_V6 = ndx_V6.difference(ndx_F3, sort=False).tolist()
F3_diff= list(set(F3_M1 + F3_L1 + F3_GVA + F3_P2 + F3_Y4 + F3_D5 + F3_H5 + F3_H6 + F3_Y6
                  + F3_Mtub2 + F3_Mtub1 +  F3_HSV1gp3 + F3_S8 + F3_V6 + F3_ImrA + F3_hCD9))

Y4_M1= ndx_M1.difference(ndx_Y4, sort=False).tolist()
Y4_L1= ndx_L1.difference(ndx_Y4, sort=False).tolist()
Y4_GVA = ndx_GVA.difference(ndx_Y4, sort=False).tolist()
Y4_P2 = ndx_P2.difference(ndx_Y4, sort=False).tolist()
Y4_F3 = ndx_F3.difference(ndx_Y4, sort=False).tolist()
Y4_D5 = ndx_D5.difference(ndx_Y4, sort=False).tolist()
Y4_H5 = ndx_H5.difference(ndx_Y4, sort=False).tolist()
Y4_H6 = ndx_H6.difference(ndx_Y4, sort=False).tolist()
Y4_Y6 = ndx_Y6.difference(ndx_Y4, sort=False).tolist()
Y4_Mtub2 = ndx_Mtub2.difference(ndx_Y4, sort=False).tolist()
Y4_ImrA = ndx_ImrA.difference(ndx_Y4, sort=False).tolist()
Y4_hCD9 = ndx_hCD9.difference(ndx_Y4, sort=False).tolist()
Y4_Mtub1 = ndx_Mtub1.difference(ndx_Y4, sort=False).tolist()
Y4_HSV1gp3 = ndx_HSV1gp3.difference(ndx_Y4, sort=False).tolist()
Y4_S8 = ndx_S8.difference(ndx_Y4, sort=False).tolist()
Y4_V6 = ndx_V6.difference(ndx_Y4, sort=False).tolist()
Y4_diff= list(set(Y4_M1 + Y4_L1 + Y4_GVA + Y4_P2 + Y4_F3 + Y4_D5 + Y4_H5 + Y4_H6 + Y4_Y6
                  + Y4_Mtub2 + Y4_Mtub1 +  Y4_HSV1gp3 + Y4_S8 + Y4_V6 + Y4_ImrA + Y4_hCD9))

D5_M1= ndx_M1.difference(ndx_D5, sort=False).tolist()
D5_L1= ndx_L1.difference(ndx_D5, sort=False).tolist()
D5_GVA = ndx_GVA.difference(ndx_D5, sort=False).tolist()
D5_P2 = ndx_P2.difference(ndx_D5, sort=False).tolist()
D5_F3 = ndx_F3.difference(ndx_D5, sort=False).tolist()
D5_Y4 = ndx_Y4.difference(ndx_D5, sort=False).tolist()
D5_H5 = ndx_H5.difference(ndx_D5, sort=False).tolist()
D5_H6 = ndx_H6.difference(ndx_D5, sort=False).tolist()
D5_Y6 = ndx_Y6.difference(ndx_D5, sort=False).tolist()
D5_Mtub2 = ndx_Mtub2.difference(ndx_D5, sort=False).tolist()
D5_ImrA = ndx_ImrA.difference(ndx_D5, sort=False).tolist()
D5_hCD9 = ndx_hCD9.difference(ndx_D5, sort=False).tolist()
D5_Mtub1 = ndx_Mtub1.difference(ndx_D5, sort=False).tolist()
D5_HSV1gp3 = ndx_HSV1gp3.difference(ndx_D5, sort=False).tolist()
D5_S8 = ndx_S8.difference(ndx_D5, sort=False).tolist()
D5_V6 = ndx_V6.difference(ndx_D5, sort=False).tolist()
D5_diff= list(set(D5_M1 + D5_L1 + D5_GVA + D5_P2 + D5_F3 + D5_Y4 + D5_H5 + D5_H6 + D5_Y6
                  + D5_Mtub2 + D5_Mtub1 +  D5_HSV1gp3 + D5_S8 + D5_V6 + D5_ImrA + D5_hCD9))

H5_M1= ndx_M1.difference(ndx_H5, sort=False).tolist()
H5_L1= ndx_L1.difference(ndx_H5, sort=False).tolist()
H5_GVA = ndx_GVA.difference(ndx_H5, sort=False).tolist()
H5_P2 = ndx_P2.difference(ndx_H5, sort=False).tolist()
H5_F3 = ndx_F3.difference(ndx_H5, sort=False).tolist()
H5_Y4 = ndx_Y4.difference(ndx_H5, sort=False).tolist()
H5_D5 = ndx_D5.difference(ndx_H5, sort=False).tolist()
H5_H6 = ndx_H6.difference(ndx_H5, sort=False).tolist()
H5_Y6 = ndx_Y6.difference(ndx_H5, sort=False).tolist()
H5_Mtub2 = ndx_Mtub2.difference(ndx_H5, sort=False).tolist()
H5_ImrA = ndx_ImrA.difference(ndx_H5, sort=False).tolist()
H5_hCD9 = ndx_hCD9.difference(ndx_H5, sort=False).tolist()
H5_Mtub1 = ndx_Mtub1.difference(ndx_H5, sort=False).tolist()
H5_HSV1gp3 = ndx_HSV1gp3.difference(ndx_H5, sort=False).tolist()
H5_S8 = ndx_S8.difference(ndx_H5, sort=False).tolist()
H5_V6 = ndx_V6.difference(ndx_H5, sort=False).tolist()
H5_diff= list(set(H5_M1 + H5_L1 + H5_GVA + H5_P2 + H5_F3 + H5_Y4 + H5_D5 + H5_H6 + H5_Y6
                  + H5_Mtub2 + H5_Mtub1 +  H5_HSV1gp3 + H5_S8 + H5_V6 + H5_ImrA + H5_hCD9))

H6_M1= ndx_M1.difference(ndx_H6, sort=False).tolist()
H6_L1= ndx_L1.difference(ndx_H6, sort=False).tolist()
H6_GVA = ndx_GVA.difference(ndx_H6, sort=False).tolist()
H6_P2 = ndx_P2.difference(ndx_H6, sort=False).tolist()
H6_F3 = ndx_F3.difference(ndx_H6, sort=False).tolist()
H6_Y4 = ndx_Y4.difference(ndx_H6, sort=False).tolist()
H6_D5 = ndx_D5.difference(ndx_H6, sort=False).tolist()
H6_H5 = ndx_H5.difference(ndx_H6, sort=False).tolist()
H6_Y6 = ndx_Y6.difference(ndx_H6, sort=False).tolist()
H6_Mtub2 = ndx_Mtub2.difference(ndx_H6, sort=False).tolist()
H6_ImrA = ndx_ImrA.difference(ndx_H6, sort=False).tolist()
H6_hCD9 = ndx_hCD9.difference(ndx_H6, sort=False).tolist()
H6_Mtub1 = ndx_Mtub1.difference(ndx_H6, sort=False).tolist()
H6_HSV1gp3 = ndx_HSV1gp3.difference(ndx_H6, sort=False).tolist()
H6_S8 = ndx_S8.difference(ndx_H6, sort=False).tolist()
H6_V6 = ndx_V6.difference(ndx_H6, sort=False).tolist()
H6_diff= list(set(H6_M1 + H6_L1 + H6_GVA + H6_P2 + H6_F3 + H6_Y4 + H6_D5 + H6_H5 + H6_Y6
                  + H6_Mtub2 + H6_Mtub1 +  H6_HSV1gp3 + H6_S8 + H6_V6 + H6_ImrA + H6_hCD9))

Y6_M1= ndx_M1.difference(ndx_Y6, sort=False).tolist()
Y6_L1= ndx_L1.difference(ndx_Y6, sort=False).tolist()
Y6_GVA = ndx_GVA.difference(ndx_Y6, sort=False).tolist()
Y6_P2 = ndx_P2.difference(ndx_Y6, sort=False).tolist()
Y6_F3 = ndx_F3.difference(ndx_Y6, sort=False).tolist()
Y6_Y4 = ndx_Y4.difference(ndx_Y6, sort=False).tolist()
Y6_D5 = ndx_D5.difference(ndx_Y6, sort=False).tolist()
Y6_H5 = ndx_H5.difference(ndx_Y6, sort=False).tolist()
Y6_H6 = ndx_H6.difference(ndx_Y6, sort=False).tolist()
Y6_Mtub2 = ndx_Mtub2.difference(ndx_Y6, sort=False).tolist()
Y6_ImrA = ndx_ImrA.difference(ndx_Y6, sort=False).tolist()
Y6_hCD9 = ndx_hCD9.difference(ndx_Y6, sort=False).tolist()
Y6_Mtub1 = ndx_Mtub1.difference(ndx_Y6, sort=False).tolist()
Y6_HSV1gp3 = ndx_HSV1gp3.difference(ndx_Y6, sort=False).tolist()
Y6_S8 = ndx_S8.difference(ndx_Y6, sort=False).tolist()
Y6_V6 = ndx_V6.difference(ndx_Y6, sort=False).tolist()
Y6_diff= list(set(Y6_M1 + Y6_L1 + Y6_GVA + Y6_P2 + Y6_F3 + Y6_Y4 + Y6_D5 + Y6_H5 + Y6_H6
                  + Y6_Mtub2 + Y6_Mtub1 +  Y6_HSV1gp3 + Y6_S8 + Y6_V6 + Y6_ImrA +Y6_hCD9))

Mtub2_M1= ndx_M1.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_L1= ndx_L1.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_GVA = ndx_GVA.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_P2 = ndx_P2.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_F3 = ndx_F3.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_Y4 = ndx_Y4.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_D5 = ndx_D5.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_H5 = ndx_H5.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_H6 = ndx_H6.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_Y6 = ndx_Y6.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_ImrA = ndx_ImrA.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_hCD9 = ndx_hCD9.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_Mtub1 = ndx_Mtub1.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_HSV1gp3 = ndx_HSV1gp3.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_S8 = ndx_S8.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_V6 = ndx_V6.difference(ndx_Mtub2, sort=False).tolist()
Mtub2_diff= list(set(Mtub2_M1 + Mtub2_L1 + Mtub2_GVA + Mtub2_P2 + Mtub2_F3 + Mtub2_Y4 + Mtub2_D5 + Mtub2_H5 + Mtub2_H6
                  + Mtub2_Y6 + Mtub2_Mtub1 +  Mtub2_HSV1gp3 + Mtub2_S8 + Mtub2_V6 + Mtub2_ImrA + Mtub2_hCD9))

Mtub1_M1= ndx_M1.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_L1= ndx_L1.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_GVA = ndx_GVA.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_P2 = ndx_P2.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_F3 = ndx_F3.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_Y4 = ndx_Y4.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_D5 = ndx_D5.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_H5 = ndx_H5.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_H6 = ndx_H6.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_Y6 = ndx_Y6.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_ImrA = ndx_ImrA.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_hCD9 = ndx_hCD9.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_Mtub2 = ndx_Mtub2.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_HSV1gp3 = ndx_HSV1gp3.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_S8 = ndx_S8.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_V6 = ndx_V6.difference(ndx_Mtub1, sort=False).tolist()
Mtub1_diff= list(set(Mtub1_M1 + Mtub1_L1 + Mtub1_GVA + Mtub1_P2 + Mtub1_F3 + Mtub1_Y4 + Mtub1_D5 + Mtub1_H5 + Mtub1_H6
                  + Mtub1_Y6 + Mtub1_Mtub2 +  Mtub1_HSV1gp3 + Mtub1_S8 + Mtub1_V6 + Mtub1_ImrA + Mtub1_hCD9))

HSV1gp3_M1= ndx_M1.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_L1= ndx_L1.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_GVA = ndx_GVA.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_P2 = ndx_P2.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_F3 = ndx_F3.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_Y4 = ndx_Y4.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_D5 = ndx_D5.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_H5 = ndx_H5.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_H6 = ndx_H6.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_Y6 = ndx_Y6.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_ImrA = ndx_ImrA.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_hCD9 = ndx_hCD9.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_Mtub2 = ndx_Mtub2.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_Mtub1 = ndx_Mtub1.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_S8 = ndx_S8.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_V6 = ndx_V6.difference(ndx_HSV1gp3, sort=False).tolist()
HSV1gp3_diff= list(set(HSV1gp3_M1 + HSV1gp3_L1 + HSV1gp3_GVA + HSV1gp3_P2 + HSV1gp3_F3 + HSV1gp3_Y4 + HSV1gp3_D5 + HSV1gp3_H5 + HSV1gp3_H6
                  + HSV1gp3_Y6 + HSV1gp3_Mtub2 +  HSV1gp3_Mtub1 + HSV1gp3_S8 + HSV1gp3_V6 + HSV1gp3_ImrA + HSV1gp3_hCD9))

S8_M1= ndx_M1.difference(ndx_S8, sort=False).tolist()
S8_L1= ndx_L1.difference(ndx_S8, sort=False).tolist()
S8_GVA = ndx_GVA.difference(ndx_S8, sort=False).tolist()
S8_P2 = ndx_P2.difference(ndx_S8, sort=False).tolist()
S8_F3 = ndx_F3.difference(ndx_S8, sort=False).tolist()
S8_Y4 = ndx_Y4.difference(ndx_S8, sort=False).tolist()
S8_D5 = ndx_D5.difference(ndx_S8, sort=False).tolist()
S8_H5 = ndx_H5.difference(ndx_S8, sort=False).tolist()
S8_H6 = ndx_H6.difference(ndx_S8, sort=False).tolist()
S8_Y6 = ndx_Y6.difference(ndx_S8, sort=False).tolist()
S8_ImrA = ndx_ImrA.difference(ndx_S8, sort=False).tolist()
S8_hCD9 = ndx_hCD9.difference(ndx_S8, sort=False).tolist()
S8_Mtub2 = ndx_Mtub2.difference(ndx_S8, sort=False).tolist()
S8_Mtub1 = ndx_Mtub1.difference(ndx_S8, sort=False).tolist()
S8_HSV1gp3 = ndx_HSV1gp3.difference(ndx_S8, sort=False).tolist()
S8_V6 = ndx_V6.difference(ndx_S8, sort=False).tolist()
S8_diff= list(set(S8_M1 + S8_L1 + S8_GVA + S8_P2 + S8_F3 + S8_Y4 + S8_D5 + S8_H5 + S8_H6
                  + S8_Y6 + S8_Mtub2 +  S8_Mtub1 + S8_HSV1gp3 + S8_V6 + S8_ImrA + S8_hCD9))

V6_M1= ndx_M1.difference(ndx_V6, sort=False).tolist()
V6_L1= ndx_L1.difference(ndx_V6, sort=False).tolist()
V6_GVA = ndx_GVA.difference(ndx_V6, sort=False).tolist()
V6_P2 = ndx_P2.difference(ndx_V6, sort=False).tolist()
V6_F3 = ndx_F3.difference(ndx_V6, sort=False).tolist()
V6_Y4 = ndx_Y4.difference(ndx_V6, sort=False).tolist()
V6_D5 = ndx_D5.difference(ndx_V6, sort=False).tolist()
V6_H5 = ndx_H5.difference(ndx_V6, sort=False).tolist()
V6_H6 = ndx_H6.difference(ndx_V6, sort=False).tolist()
V6_Y6 = ndx_Y6.difference(ndx_V6, sort=False).tolist()
V6_ImrA = ndx_ImrA.difference(ndx_V6, sort=False).tolist()
V6_hCD9 = ndx_hCD9.difference(ndx_V6, sort=False).tolist()
V6_Mtub2 = ndx_Mtub2.difference(ndx_V6, sort=False).tolist()
V6_Mtub1 = ndx_Mtub1.difference(ndx_V6, sort=False).tolist()
V6_HSV1gp3 = ndx_HSV1gp3.difference(ndx_V6, sort=False).tolist()
V6_S8 = ndx_S8.difference(ndx_V6, sort=False).tolist()
V6_diff= list(set(V6_M1 + V6_L1 + V6_GVA + V6_P2 + V6_F3 + V6_Y4 + V6_D5 + V6_H5 + V6_H6
                  + V6_Y6 + V6_Mtub2 +  V6_Mtub1 + V6_HSV1gp3 + V6_S8 + V6_ImrA + V6_hCD9))

ImrA_M1= ndx_M1.difference(ndx_ImrA, sort=False).tolist()
ImrA_L1= ndx_L1.difference(ndx_ImrA, sort=False).tolist()
ImrA_GVA = ndx_GVA.difference(ndx_ImrA, sort=False).tolist()
ImrA_P2 = ndx_P2.difference(ndx_ImrA, sort=False).tolist()
ImrA_F3 = ndx_F3.difference(ndx_ImrA, sort=False).tolist()
ImrA_Y4 = ndx_Y4.difference(ndx_ImrA, sort=False).tolist()
ImrA_D5 = ndx_D5.difference(ndx_ImrA, sort=False).tolist()
ImrA_H5 = ndx_H5.difference(ndx_ImrA, sort=False).tolist()
ImrA_H6 = ndx_H6.difference(ndx_ImrA, sort=False).tolist()
ImrA_Y6 = ndx_Y6.difference(ndx_ImrA, sort=False).tolist()
ImrA_V6 = ndx_V6.difference(ndx_ImrA, sort=False).tolist()
ImrA_hCD9 = ndx_hCD9.difference(ndx_ImrA, sort=False).tolist()
ImrA_Mtub2 = ndx_Mtub2.difference(ndx_ImrA, sort=False).tolist()
ImrA_Mtub1 = ndx_Mtub1.difference(ndx_ImrA, sort=False).tolist()
ImrA_HSV1gp3 = ndx_HSV1gp3.difference(ndx_ImrA, sort=False).tolist()
ImrA_S8 = ndx_S8.difference(ndx_ImrA, sort=False).tolist()
ImrA_diff= list(set(ImrA_M1 + ImrA_L1 + ImrA_GVA + ImrA_P2 + ImrA_F3 + ImrA_Y4 + ImrA_D5 + ImrA_H5 + ImrA_H6
                  + ImrA_Y6 + ImrA_Mtub2 +  ImrA_Mtub1 + ImrA_HSV1gp3 + ImrA_S8 + ImrA_V6 + ImrA_hCD9))


hCD9_M1= ndx_M1.difference(ndx_hCD9, sort=False).tolist()
hCD9_L1= ndx_L1.difference(ndx_hCD9, sort=False).tolist()
hCD9_GVA = ndx_GVA.difference(ndx_hCD9, sort=False).tolist()
hCD9_P2 = ndx_P2.difference(ndx_hCD9, sort=False).tolist()
hCD9_F3 = ndx_F3.difference(ndx_hCD9, sort=False).tolist()
hCD9_Y4 = ndx_Y4.difference(ndx_hCD9, sort=False).tolist()
hCD9_D5 = ndx_D5.difference(ndx_hCD9, sort=False).tolist()
hCD9_H5 = ndx_H5.difference(ndx_hCD9, sort=False).tolist()
hCD9_H6 = ndx_H6.difference(ndx_hCD9, sort=False).tolist()
hCD9_Y6 = ndx_Y6.difference(ndx_hCD9, sort=False).tolist()
hCD9_V6 = ndx_V6.difference(ndx_hCD9, sort=False).tolist()
hCD9_ImrA = ndx_ImrA.difference(ndx_hCD9, sort=False).tolist()
hCD9_Mtub2 = ndx_Mtub2.difference(ndx_hCD9, sort=False).tolist()
hCD9_Mtub1 = ndx_Mtub1.difference(ndx_hCD9, sort=False).tolist()
hCD9_HSV1gp3 = ndx_HSV1gp3.difference(ndx_hCD9, sort=False).tolist()
hCD9_S8 = ndx_S8.difference(ndx_hCD9, sort=False).tolist()
hCD9_diff= list(set(hCD9_M1 + hCD9_L1 + hCD9_GVA + hCD9_P2 + hCD9_F3 + hCD9_Y4 + hCD9_D5 + hCD9_H5 + hCD9_H6
                  + hCD9_Y6 + hCD9_Mtub2 +  hCD9_Mtub1 + hCD9_HSV1gp3 + hCD9_S8 + hCD9_V6 + hCD9_ImrA))

###Inter-Mutant Total Common Index Addition Loops###

#L1#



print(len(L1))
print(len(L1_diff))
for i in range(len(L1_diff)):
    add_zeros = np.zeros(shape=(1, len(L1.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [L1_diff[i]])
    L1 = L1.append(zero)

L1 = L1.sort_index(axis=0, ascending=True)
print(len(L1))


#MART1#
print(len(M1_diff))
print(len(MART1))
for i in range(len(M1_diff)):
    add_zeros = np.zeros(shape=(1, len(MART1.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [M1_diff[i]])
    MART1 = MART1.append(zero)

MART1 = MART1.sort_index(axis=0, ascending=True)
print(len(MART1))

#GVA#
print(len(GVA_diff))
print(len(GVA))
for i in range(len(GVA_diff)):
    add_zeros = np.zeros(shape=(1, len(GVA.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [GVA_diff[i]])
    GVA = GVA.append(zero)
 
GVA = GVA.sort_index(axis=0, ascending=True)
print(len(GVA))

#P2#
print(len(P2_diff))
print(len(P2))
for i in range(len(P2_diff)):
    add_zeros = np.zeros(shape=(1, len(P2.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [P2_diff[i]])
    P2 = P2.append(zero)
 
P2 = P2.sort_index(axis=0, ascending=True)
print(len(P2))

#F3#
print(len(F3_diff))
print(len(F3))
for i in range(len(F3_diff)):
    add_zeros = np.zeros(shape=(1, len(F3.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [F3_diff[i]])
    F3 = F3.append(zero)
 
F3 = F3.sort_index(axis=0, ascending=True)
print(len(F3))

#Y4#
print(len(Y4_diff))
print(len(Y4))
for i in range(len(Y4_diff)):
    add_zeros = np.zeros(shape=(1, len(Y4.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [Y4_diff[i]])
    Y4 = Y4.append(zero)
 
Y4 = Y4.sort_index(axis=0, ascending=True)
print(len(Y4))

#D5#
print(len(D5_diff))
print(len(D5))
for i in range(len(D5_diff)):
    add_zeros = np.zeros(shape=(1, len(D5.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [D5_diff[i]])
    D5 = D5.append(zero)
 
D5 = D5.sort_index(axis=0, ascending=True)
print(len(D5))

#H5#
print(len(H5_diff))
print(len(H5))
for i in range(len(H5_diff)):
    add_zeros = np.zeros(shape=(1, len(H5.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [H5_diff[i]])
    H5 = H5.append(zero)
 
H5 = H5.sort_index(axis=0, ascending=True)
print(len(H5))

#H6#
print(len(H6_diff))
print(len(H6))
for i in range(len(H6_diff)):
    add_zeros = np.zeros(shape=(1, len(H6.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [H6_diff[i]])
    H6 = H6.append(zero)
 
H6 = H6.sort_index(axis=0, ascending=True)
print(len(H6))

#Y6#
print(len(Y6_diff))
print(len(Y6))
for i in range(len(Y6_diff)):
    add_zeros = np.zeros(shape=(1, len(Y6.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [Y6_diff[i]])
    Y6 = Y6.append(zero)
 
Y6 = Y6.sort_index(axis=0, ascending=True)
print(len(Y6))

#Mtub2#
print(len(Mtub2_diff))
print(len(Mtub2))
for i in range(len(Mtub2_diff)):
    add_zeros = np.zeros(shape=(1, len(Mtub2.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [Mtub2_diff[i]])
    Mtub2 = Mtub2.append(zero)
 
Mtub2 = Mtub2.sort_index(axis=0, ascending=True)
print(len(Mtub2))

#Mtub1#
print(len(Mtub1_diff))
print(len(Mtub1))
for i in range(len(Mtub1_diff)):
    add_zeros = np.zeros(shape=(1, len(Mtub1.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [Mtub1_diff[i]])
    Mtub1 = Mtub1.append(zero)
 
Mtub1 = Mtub1.sort_index(axis=0, ascending=True)
print(len(Mtub1))

#HSV1gp3#
print(len(HSV1gp3_diff))
print(len(HSV1gp3))
for i in range(len(HSV1gp3_diff)):
    add_zeros = np.zeros(shape=(1, len(HSV1gp3.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [HSV1gp3_diff[i]])
    HSV1gp3 = HSV1gp3.append(zero)
 
HSV1gp3 = HSV1gp3.sort_index(axis=0, ascending=True)
print(len(HSV1gp3))

#S8#
print(len(S8_diff))
print(len(S8))
for i in range(len(S8_diff)):
    add_zeros = np.zeros(shape=(1, len(S8.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [S8_diff[i]])
    S8 = S8.append(zero)
 
S8 = S8.sort_index(axis=0, ascending=True)
print(len(S8))

#V6#
print(len(V6_diff))
print(len(V6))
for i in range(len(V6_diff)):
    add_zeros = np.zeros(shape=(1, len(V6.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [V6_diff[i]])
    V6 = V6.append(zero)
 
V6 = V6.sort_index(axis=0, ascending=True)
print(len(V6))

#ImrA#
print(len(ImrA_diff))
print(len(ImrA))
for i in range(len(ImrA_diff)):
    add_zeros = np.zeros(shape=(1, len(ImrA.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [ImrA_diff[i]])
    ImrA = ImrA.append(zero)
 
ImrA = ImrA.sort_index(axis=0, ascending=True)
print(len(ImrA))

#hCD9#
print(len(hCD9_diff))
print(len(hCD9))
for i in range(len(hCD9_diff)):
    add_zeros = np.zeros(shape=(1, len(hCD9.columns)-1))
    zero = pd.DataFrame(add_zeros, index = [hCD9_diff[i]])
    hCD9 = hCD9.append(zero)
 
hCD9 = hCD9.sort_index(axis=0, ascending=True)
print(len(hCD9))

#FillNaN

L1 = L1.fillna(0)
MART1 = MART1.fillna(0)
GVA = GVA.fillna(0)
P2 = P2.fillna(0)
F3 = F3.fillna(0)
Y4 = Y4.fillna(0)
D5 = D5.fillna(0)
H5 = H5.fillna(0)
H6 = H6.fillna(0)
Y6 = Y6.fillna(0)
Mtub2 = Mtub2.fillna(0)
Mtub1 = Mtub1.fillna(0)
HSV1gp3 = HSV1gp3.fillna(0)
S8 = S8.fillna(0)
V6 = V6.fillna(0)
ImrA= ImrA.fillna(0)
hCD9 = hCD9.fillna(0)


#Save Master Maps

L1_MHC = L1[L1.index.str.contains('MHC')]
L1.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/hmapx_pmhc_master.csv', header=None)
L1_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/hmapx_mhc_master.csv', header=None)

MART1_MHC = MART1[MART1.index.str.contains('MHC')]
MART1.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/hmapx_pmhc_master.csv', header=None)
MART1_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/hmapx_mhc_master.csv', header=None)

GVA_MHC = GVA[GVA.index.str.contains('MHC')]
GVA.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/hmapx_pmhc_master.csv', header=None)
GVA_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/hmapx_mhc_master.csv', header=None)

P2_MHC = P2[P2.index.str.contains('MHC')]
P2.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/2P/hmapx_pmhc_master.csv', header=None)
P2_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/2P/hmapx_mhc_master.csv', header=None)

F3_MHC = F3[F3.index.str.contains('MHC')]
F3.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/3F/hmapx_pmhc_master.csv', header=None)
F3_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/3F/hmapx_mhc_master.csv', header=None)

D5_MHC = D5[D5.index.str.contains('MHC')]
D5.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5D/hmapx_pmhc_master.csv', header=None)
D5_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5D/hmapx_mhc_master.csv', header=None)

Y4_MHC = Y4[Y4.index.str.contains('MHC')]
Y4.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/4Y/hmapx_pmhc_master.csv', header=None)
Y4_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/4Y/hmapx_mhc_master.csv', header=None)

H5_MHC = H5[H5.index.str.contains('MHC')]
H5.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5H/hmapx_pmhc_master.csv', header=None)
H5_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5H/hmapx_mhc_master.csv', header=None)

H6_MHC = H6[H6.index.str.contains('MHC')]
H6.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6H/hmapx_pmhc_master.csv', header=None)
H6_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6H/hmapx_mhc_master.csv', header=None)

Y6_MHC = Y6[Y6.index.str.contains('MHC')]
Y6.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6Y/hmapx_pmhc_master.csv', header=None)
Y6_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6Y/hmapx_mhc_master.csv', header=None)

Mtub2_MHC = Mtub2[Mtub2.index.str.contains('MHC')]
Mtub2.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub2/hmapx_pmhc_master.csv', header=None)
Mtub2_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub2/hmapx_mhc_master.csv', header=None)

Mtub1_MHC = Mtub1[Mtub1.index.str.contains('MHC')]
Mtub1.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub1/hmapx_pmhc_master.csv', header=None)
Mtub1_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub1/hmapx_mhc_master.csv', header=None)

HSV1gp3_MHC = HSV1gp3[HSV1gp3.index.str.contains('MHC')]
HSV1gp3.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/HSV1gp3/hmapx_pmhc_master.csv', header=None)
HSV1gp3_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/HSV1gp3/hmapx_mhc_master.csv', header=None)

S8_MHC = S8[S8.index.str.contains('MHC')]
S8.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/hmapx_pmhc_master.csv', header=None)
S8_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/hmapx_mhc_master.csv', header=None)

V6_MHC = V6[V6.index.str.contains('MHC')]
V6.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/6V/hmapx_pmhc_master.csv', header=None)
V6_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/6V/hmapx_mhc_master.csv', header=None)

ImrA_MHC = ImrA[ImrA.index.str.contains('MHC')]
ImrA.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/ImrA/hmapx_pmhc_master.csv', header=None)
ImrA_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/ImrA/hmapx_mhc_master.csv', header=None)

hCD9_MHC = hCD9[hCD9.index.str.contains('MHC')]
hCD9.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/hmapx_pmhc_master.csv', header=None)
hCD9_MHC.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/hmapx_mhc_master.csv', header=None)
