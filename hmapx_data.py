#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 18:42:34 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch


#Interfacial Residues
resi = ['CDR1a_1','CDR1a_2','CDR1a_3','CDR1a_4','CDR1a_5','CDR1a_6','CDR1a_7','CDR1a_8','CDR1a_9',
        'CDR2a_1','CDR2a_2','CDR2a_3','CDR2a_4','CDR2a_5','CDR2a_6',
        'CDR3a_1','CDR3a_2','CDR3a_3','CDR3a_4','CDR3a_5','CDR3a_6','CDR3a_7','CDR3a_8','CDR3a_9','CDR3a_10','CDR3a_11',
        'CDR1b_1','CDR1b_2','CDR1b_3','CDR1b_4','CDR1b_5','CDR1b_6','CDR1b_7',
        'CDR2b_1','CDR2b_2','CDR2b_3','CDR2b_4','CDR2b_5','CDR2b_6','CDR2b_7','CDR2b_8',
        'CDR3b_1','CDR3b_2','CDR3b_3','CDR3b_4','CDR3b_5','CDR3b_6','CDR3b_7','CDR3b_8','CDR3b_9','CDR3b_10','CDR3b_11','CDR3b_12',
        'MHCa_1','MHCa_2','MHCa_3','MHCa_4','MHCa_5','MHCa_6','MHCa_7','MHCa_8','MHCa_9','MHCa_10','MHCa_11','MHCa_12',
        'MHCa_13','MHCa_14','MHCa_15','MHCa_16','MHCa_17','MHCa_18','MHCa_19','MHCa_20','MHCa_21','MHCa_22','MHCa_23','MHCa_24',
        'MHCa_25','MHCa_26','MHCa_27','MHCa_28','MHCa_29','MHCa_30','MHCa_31','MHCa_32','MHCa_33','MHCa_34','MHCa_35','MHCa_36',
        'MHCa_37','MHCa_38','MHCa_39','MHCa_40','MHCa_41','MHCa_42',
        'MHCb_1','MHCb_2','MHCb_3','MHCb_4','MHCb_5','MHCb_6','MHCb_7','MHCb_8','MHCb_9','MHCb_10','MHCb_11','MHCb_12',
        'MHCb_13','MHCb_14','MHCb_15','MHCb_16','MHCb_17','MHCb_18','MHCb_19','MHCb_20','MHCb_21','MHCb_22','MHCb_23','MHCb_24',
        'MHCb_25','MHCb_26','MHCb_27','MHCb_28','MHCb_29','MHCb_30','MHCb_31','MHCb_32','MHCb_33','MHCb_34','MHCb_35','MHCb_36',
        'pep_1','pep_2','pep_3','pep_4','pep_5','pep_6','pep_7','pep_8','pep_9']

#Transform Data for RNN, add resi
resi_L1 = ['Y','S','D','R','G','S','Q','S','F',
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
        'L','A','G','I','G','I','L','T','V']

L1 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/hmapx_pmhc_master.csv', index_col = 0, header=None)
L1 = L1.stack(dropna=False)
print(L1)
print(L1.shape)

resi_M1 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','G','I','L','T','V']

M1 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/hmapx_pmhc_master.csv', index_col = 0, header=None)
M1 = M1.stack(dropna=False)
print(M1)
print(M1.shape)

resi_GVA = ['Y','S','D','R','G','S','Q','S','F',
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
        'G','A','G','I','G','V','L','T','A']

GVA = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/hmapx_pmhc_master.csv', index_col = 0, header=None)
GVA = GVA.stack(dropna=False)
print(GVA)
print(GVA.shape)

resi_P2 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','P','G','I','G','I','L','T','V']

P2 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/2P/hmapx_pmhc_master.csv', index_col = 0, header=None)
P2 = P2.stack(dropna=False)
print(P2)
print(P2.shape)

resi_F3 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','F','I','G','I','L','T','V']

F3 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/3F/hmapx_pmhc_master.csv', index_col = 0, header=None)
F3 = F3.stack(dropna=False)
print(F3)
print(F3.shape)

resi_Y4 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','Y','G','I','L','T','V']

Y4 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/4Y/hmapx_pmhc_master.csv', index_col = 0, header=None)
Y4 = Y4.stack(dropna=False)
print(Y4)
print(Y4.shape)

resi_D5 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','D','I','L','T','V']

D5 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5D/hmapx_pmhc_master.csv', index_col = 0, header=None)
D5 = D5.stack(dropna=False)
print(D5)
print(D5.shape)

resi_H5 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','H','I','L','T','V']

H5 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/5H/hmapx_pmhc_master.csv', index_col = 0, header=None)
H5 = H5.stack(dropna=False)
print(H5)
print(H5.shape)

resi_H6 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','G','H','L','T','V']

H6 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6H/hmapx_pmhc_master.csv', index_col = 0, header=None)
H6 = H6.stack(dropna=False)
print(H6)
print(H6.shape)

resi_Y6 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','G','Y','L','T','V']

Y6 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/6Y/hmapx_pmhc_master.csv', index_col = 0, header=None)
Y6 = Y6.stack(dropna=False)
print(Y6)
print(Y6.shape)

resi_Mtub2 = ['Y','S','D','R','G','S','Q','S','F',
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
        'I','A','G','P','G','T','I','T','L']

Mtub2 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub2/hmapx_pmhc_master.csv', index_col = 0, header=None)
Mtub2 = Mtub2.stack(dropna=False)
print(Mtub2)
print(Mtub2.shape)

resi_Mtub1 = ['Y','S','D','R','G','S','Q','S','F',
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
        'L','G','G','L','G','L','F','F','A']

Mtub1 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/Mtub1/hmapx_pmhc_master.csv', index_col = 0, header=None)
Mtub1 = Mtub1.stack(dropna=False)
print(Mtub1)
print(Mtub1.shape)

resi_HSV1gp3 = ['Y','S','D','R','G','S','Q','S','F',
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
        'I','A','G','I','G','I','L','A','I']

HSV1gp3 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/HSV1gp3/hmapx_pmhc_master.csv', index_col = 0, header=None)
HSV1gp3 = HSV1gp3.stack(dropna=False)
print(HSV1gp3)
print(HSV1gp3.shape)

resi_S8 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','G','I','L','S','V']

S8 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/hmapx_pmhc_master.csv', index_col = 0, header=None)
S8 = S8.stack(dropna=False)
print(S8)
print(S8.shape)

resi_V6 = ['Y','S','D','R','G','S','Q','S','F',
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
        'A','A','G','I','G','V','L','T','V']

V6 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/6V/hmapx_pmhc_master.csv', index_col = 0, header=None)
V6 = V6.stack(dropna=False)
print(V6)
print(V6.shape)

resi_ImrA = ['Y','S','D','R','G','S','Q','S','F',
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
        'L','A','G','I','G','L','I','A','A']

ImrA = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/ImrA/hmapx_pmhc_master.csv', index_col = 0, header=None)
ImrA = ImrA.stack(dropna=False)
print(ImrA)
print(ImrA.shape)

resi_hCD9 = ['Y','S','D','R','G','S','Q','S','F',
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

hCD9 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/hmapx_pmhc_master.csv', index_col = 0, header=None)
hCD9 = hCD9.stack(dropna=False)
print(hCD9)
print(hCD9.shape)

#Compile DataTransform for RNN

datar_labels = {'L1':resi_L1,'MART1':resi_M1,'GVA':resi_GVA
                ,'Mtub2':resi_Mtub2, 'Mtub1':resi_Mtub1,'HSV1gp3':resi_HSV1gp3, '8S':resi_S8,'6V':resi_V6, 'ImrA':resi_ImrA,'hCD9':resi_hCD9}
#,'2P':resi_P2,'3F':resi_F3,'4Y':resi_Y4,'5D':resi_D5,'5H':resi_H5,'6H':resi_H6,'6Y':resi_Y6
datar = pd.DataFrame(datar_labels)
print(datar)
datar = datar.replace(['A','R','N','D','C','Q','E','G','H','I','L',
                         'K','M','F','P','S','T','W','Y','V']
                       ,['1','2','3','4','5','6','7','8','9','10',
                       '11','12','13','14','15','16','17','18','19','20'])


data_labels = {'L1':L1,'MART1':M1,'GVA':GVA,'Mtub2':Mtub2, 'Mtub1':Mtub1,'HSV1gp3':HSV1gp3,'8S':S8,'6V':V6, 'ImrA':ImrA,'hCD9':hCD9}

               #'Mtub2':Mtub2, 'Mtub1':Mtub1,'HSV1gp3':HSV1gp3,'8S':S8,'6V':V6, 'ImrA':ImrA,'hCD9':hCD9
               #,'2P':P2,'3F':F3,'4Y':Y4,'5D':D5,'5H':H5,'6H':H6,'6Y':Y6
data = pd.DataFrame(data_labels)
print(data)

datar = datar.append(data, ignore_index=True)
datar = datar.T

print(datar)


datar.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/hmapx_pmhc_master_bio.csv', header=None)



