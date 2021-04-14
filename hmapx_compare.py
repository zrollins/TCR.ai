#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 09:48:48 2021

@author: zrollins
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch


#Stacked Data
data = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/hmapx_pmhc_master.csv', index_col = 0, header=None)
data = data.iloc[:3,140:].T.reset_index(drop=True)
data.index = data.index +1
data = data.T
L1 = data.iloc[:1,:].T
M1 = data.iloc[1:2,:].T
GVA = data.iloc[2:3,:].T

#Stacked Predictions
prd = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/hmapx_pmhc_mstr_prd_2x128_rndm.csv', index_col = 0, header=None)
L1_prd = prd.iloc[:1,:].T
M1_prd = prd.iloc[1:2,:].T
GVA_prd = prd.iloc[2:3,:].T

#Stacked Indices
L1_map = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/hmapx_pmhc_master.csv', index_col = 0, header=None)
L1_map = L1_map.stack(dropna=False)
M1_map = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/hmapx_pmhc_master.csv', index_col = 0, header=None)
M1_map = M1_map.stack(dropna=False)
GVA_map = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/hmapx_pmhc_master.csv', index_col = 0, header=None)
GVA_map = GVA_map.stack(dropna=False)

#Set Stack Indices on Prediction and Data
L1_prd.set_index(L1_map.index, inplace=True)
M1_prd.set_index(M1_map.index, inplace=True)
GVA_prd.set_index(GVA_map.index, inplace=True)

L1.set_index(L1_map.index, inplace=True)
M1.set_index(M1_map.index, inplace=True)
GVA.set_index(GVA_map.index, inplace=True)

#Unstack
L1_prd = L1_prd.unstack()
M1_prd = M1_prd.unstack()
GVA_prd = GVA_prd.unstack()

L1 = L1.unstack()
M1 = M1.unstack()
GVA = GVA.unstack()

# High Occupancy Filter
thresh = 0.05
L1['AVG Fractional Occupancy (%)'] = L1.mean(axis=1).round(decimals=3)
L1['AVG Fractional Occupancy (%)'] = 100*L1['AVG Fractional Occupancy (%)']
L1 = L1[L1['AVG Fractional Occupancy (%)']> thresh]
print(L1.shape)

M1['AVG Fractional Occupancy (%)'] = M1.mean(axis=1).round(decimals=3)
M1['AVG Fractional Occupancy (%)'] = 100*M1['AVG Fractional Occupancy (%)']
M1 = M1[M1['AVG Fractional Occupancy (%)']> thresh]
print(M1.shape)

GVA['AVG Fractional Occupancy (%)'] = GVA.mean(axis=1).round(decimals=3)
GVA['AVG Fractional Occupancy (%)'] = 100*GVA['AVG Fractional Occupancy (%)']
GVA = GVA[GVA['AVG Fractional Occupancy (%)']> thresh]
print(GVA.shape)

#Match Data and Prediction Indices
L1_prd = L1_prd.loc[L1.index,:]
L1_prd['AVG Fractional Occupancy (%)'] = L1_prd.mean(axis=1).round(decimals=3)
L1_prd['AVG Fractional Occupancy (%)'] = 100*L1_prd['AVG Fractional Occupancy (%)']
M1_prd = M1_prd.loc[M1.index,:]
M1_prd['AVG Fractional Occupancy (%)'] = M1_prd.mean(axis=1).round(decimals=3)
M1_prd['AVG Fractional Occupancy (%)'] = 100*M1_prd['AVG Fractional Occupancy (%)']
GVA_prd = GVA_prd.loc[GVA.index,:]
GVA_prd['AVG Fractional Occupancy (%)'] = GVA_prd.mean(axis=1).round(decimals=3)
GVA_prd['AVG Fractional Occupancy (%)'] = 100*GVA_prd['AVG Fractional Occupancy (%)']

#Drop Occupancy
L1_ocpy = pd.DataFrame()
L1_ocpy_prd = pd.DataFrame()
L1_ocpy['AVG Fractional Occupancy (%)'] = L1['AVG Fractional Occupancy (%)']
L1 = L1.drop('AVG Fractional Occupancy (%)', axis=1)
L1_ocpy_prd['AVG Fractional Occupancy (%)'] = L1_prd['AVG Fractional Occupancy (%)']
L1_prd = L1_prd.drop('AVG Fractional Occupancy (%)', axis=1)

M1_ocpy = pd.DataFrame()
M1_ocpy_prd = pd.DataFrame()
M1_ocpy['AVG Fractional Occupancy (%)'] = M1['AVG Fractional Occupancy (%)']
M1 = M1.drop('AVG Fractional Occupancy (%)', axis=1)
M1_ocpy_prd['AVG Fractional Occupancy (%)'] = M1_prd['AVG Fractional Occupancy (%)']
M1_prd = M1_prd.drop('AVG Fractional Occupancy (%)', axis=1)

GVA_ocpy = pd.DataFrame()
GVA_ocpy_prd = pd.DataFrame()
GVA_ocpy['AVG Fractional Occupancy (%)'] = GVA['AVG Fractional Occupancy (%)']
GVA = GVA.drop('AVG Fractional Occupancy (%)', axis=1)
GVA_ocpy_prd['AVG Fractional Occupancy (%)'] = GVA_prd['AVG Fractional Occupancy (%)']
GVA_prd = GVA_prd.drop('AVG Fractional Occupancy (%)', axis=1)

#Separate Peptide and MHC H-Bonds
L1_MHC = L1[L1.index.str.contains('MHC')]
L1_pep = L1[L1.index.str.contains('L1')]
L1_MHC_prd = L1_prd[L1_prd.index.str.contains('MHC')]
L1_pep_prd = L1_prd[L1_prd.index.str.contains('L1')]

M1_MHC = M1[M1.index.str.contains('MHC')]
M1_pep = M1[M1.index.str.contains('MART1')]
M1_MHC_prd = M1_prd[M1_prd.index.str.contains('MHC')]
M1_pep_prd = M1_prd[M1_prd.index.str.contains('MART1')]

GVA_MHC = GVA[GVA.index.str.contains('MHC')]
GVA_pep = GVA[GVA.index.str.contains('GVA')]
GVA_MHC_prd = GVA_prd[GVA_prd.index.str.contains('MHC')]
GVA_pep_prd = GVA_prd[GVA_prd.index.str.contains('GVA')]

# Plot H-Bond Maps (Data vs. RNN Prediction)

FS=7
bins=np.linspace(5.884,9.000,200).tolist()
xticks=bins[::20]
xticks=np.around(xticks,2)

#L1
MHC_rows, MHC_columns = L1_MHC.shape
pep_rows, pep_columns = L1_pep.shape
fig, ax = plt.subplots(2,sharex=True,sharey=False,figsize=(20, 10),gridspec_kw={'height_ratios':[MHC_rows,4]})
plot=ax[0].imshow(L1_MHC,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[0].set_yticks(range(len(L1_MHC.index)))
ax[0].set_yticklabels(L1_MHC.index, fontsize=FS)
ax[0].set_title('High Occupancy H-Bond Heat Map (L1)',fontsize=20)

ax[1].imshow(L1_pep,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[1].set_yticks(range(len(L1_pep.index)))
ax[1].set_yticklabels(L1_pep.index, fontsize=FS)
ax[1].set_xlabel('Reaction Coordinate, ξ (nm)', fontsize=15)
ax[1].set_xticks(range(0,len(L1_pep.columns.tolist()),20))
ax[1].set_xticklabels(xticks)

fig.colorbar(plot, ax=ax[0:], location='right')
fig.text(-0.03,0.5, 'Hydrogen Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/HBondx.png', bbox_inches='tight', dpi=300)

MHC_rows, MHC_columns = L1_MHC_prd.shape
pep_rows, pep_columns = L1_pep_prd.shape
fig, ax = plt.subplots(2,sharex=True,sharey=False,figsize=(20, 10),gridspec_kw={'height_ratios':[MHC_rows,4]})
plot=ax[0].imshow(L1_MHC_prd,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[0].set_yticks(range(len(L1_MHC_prd.index)))
ax[0].set_yticklabels(L1_MHC_prd.index, fontsize=FS)
ax[0].set_title('High Occupancy H-Bond Heat Map RNN Prediction (L1)',fontsize=20)

ax[1].imshow(L1_pep_prd,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[1].set_yticks(range(len(L1_pep_prd.index)))
ax[1].set_yticklabels(L1_pep_prd.index, fontsize=FS)
ax[1].set_xlabel('Reaction Coordinate, ξ (nm)', fontsize=15)
ax[1].set_xticks(range(0,len(L1_pep_prd.columns.tolist()),20))
ax[1].set_xticklabels(xticks)

fig.colorbar(plot, ax=ax[0:], location='right')
fig.text(-0.03,0.5, 'Hydrogen Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/HBondx_RNN_2x128_rndm.png', bbox_inches='tight', dpi=300)

#M1
MHC_rows, MHC_columns = M1_MHC.shape
pep_rows, pep_columns = M1_pep.shape
fig, ax = plt.subplots(2,sharex=True,sharey=False,figsize=(20, 10),gridspec_kw={'height_ratios':[MHC_rows,4]})
plot=ax[0].imshow(M1_MHC,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[0].set_yticks(range(len(M1_MHC.index)))
ax[0].set_yticklabels(M1_MHC.index, fontsize=FS)
ax[0].set_title('High Occupancy H-Bond Heat Map (M1)',fontsize=20)

ax[1].imshow(M1_pep,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[1].set_yticks(range(len(M1_pep.index)))
ax[1].set_yticklabels(M1_pep.index, fontsize=FS)
ax[1].set_xlabel('Reaction Coordinate, ξ (nm)', fontsize=15)
ax[1].set_xticks(range(0,len(M1_pep.columns.tolist()),20))
ax[1].set_xticklabels(xticks)

fig.colorbar(plot, ax=ax[0:], location='right')
fig.text(-0.03,0.5, 'Hydrogen Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/HBondx.png', bbox_inches='tight', dpi=300)

MHC_rows, MHC_columns = M1_MHC_prd.shape
pep_rows, pep_columns = M1_pep_prd.shape
fig, ax = plt.subplots(2,sharex=True,sharey=False,figsize=(20, 10),gridspec_kw={'height_ratios':[MHC_rows,4]})
plot=ax[0].imshow(M1_MHC_prd,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[0].set_yticks(range(len(M1_MHC_prd.index)))
ax[0].set_yticklabels(M1_MHC_prd.index, fontsize=FS)
ax[0].set_title('High Occupancy H-Bond Heat Map RNN Prediction (M1)',fontsize=20)

ax[1].imshow(M1_pep_prd,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[1].set_yticks(range(len(M1_pep_prd.index)))
ax[1].set_yticklabels(M1_pep_prd.index, fontsize=FS)
ax[1].set_xlabel('Reaction Coordinate, ξ (nm)', fontsize=15)
ax[1].set_xticks(range(0,len(M1_pep_prd.columns.tolist()),20))
ax[1].set_xticklabels(xticks)

fig.colorbar(plot, ax=ax[0:], location='right')
fig.text(-0.03,0.5, 'Hydrogen Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/HBondx_RNN_2x128_rndm.png', bbox_inches='tight', dpi=300)

#GVA
MHC_rows, MHC_columns = GVA_MHC.shape
pep_rows, pep_columns = GVA_pep.shape
fig, ax = plt.subplots(2,sharex=True,sharey=False,figsize=(20, 10),gridspec_kw={'height_ratios':[MHC_rows,4]})
plot=ax[0].imshow(GVA_MHC,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[0].set_yticks(range(len(GVA_MHC.index)))
ax[0].set_yticklabels(GVA_MHC.index, fontsize=FS)
ax[0].set_title('High Occupancy H-Bond Heat Map (GVA)',fontsize=20)

ax[1].imshow(GVA_pep,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[1].set_yticks(range(len(GVA_pep.index)))
ax[1].set_yticklabels(GVA_pep.index, fontsize=FS)
ax[1].set_xlabel('Reaction Coordinate, ξ (nm)', fontsize=15)
ax[1].set_xticks(range(0,len(GVA_pep.columns.tolist()),20))
ax[1].set_xticklabels(xticks)

fig.colorbar(plot, ax=ax[0:], location='right')
fig.text(-0.03,0.5, 'Hydrogen Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/HBondx.png', bbox_inches='tight', dpi=300)

MHC_rows, MHC_columns = GVA_MHC_prd.shape
pep_rows, pep_columns = GVA_pep_prd.shape
fig, ax = plt.subplots(2,sharex=True,sharey=False,figsize=(20, 10),gridspec_kw={'height_ratios':[MHC_rows,4]})
plot=ax[0].imshow(GVA_MHC_prd,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[0].set_yticks(range(len(GVA_MHC_prd.index)))
ax[0].set_yticklabels(GVA_MHC_prd.index, fontsize=FS)
ax[0].set_title('High Occupancy H-Bond Heat Map RNN Prediction (GVA)',fontsize=20)

ax[1].imshow(GVA_pep_prd,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
ax[1].set_yticks(range(len(GVA_pep_prd.index)))
ax[1].set_yticklabels(GVA_pep_prd.index, fontsize=FS)
ax[1].set_xlabel('Reaction Coordinate, ξ (nm)', fontsize=15)
ax[1].set_xticks(range(0,len(GVA_pep_prd.columns.tolist()),20))
ax[1].set_xticklabels(xticks)

fig.colorbar(plot, ax=ax[0:], location='right')
fig.text(-0.03,0.5, 'Hydrogen Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/HBondx_RNN_2x128_rndm.png', bbox_inches='tight', dpi=300)



