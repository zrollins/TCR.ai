#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 09:16:46 2021

@author: zrollins
"""

import numpy as np
from numpy import mean
from numpy import std
from numpy import absolute
import pandas as pd
import matplotlib.pyplot as plt
from numpy import arange






pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 1000)

q = pd.read_excel(r'/Users/zrollins/Documents/Documents/TCR_ML/InformationCriteria.xlsx', sheet_name='qfeatures', index_col = 0, header=None)


q_ic = pd.DataFrame()
q_ic1 = pd.DataFrame()
q_ic3 = pd.DataFrame()
q_ic5= pd.DataFrame()
q_ic7 = pd.DataFrame()
q_ic['Model']=q.iloc[1:9,0].reset_index(drop=True)
q_ic1['AIC'] = q.iloc[1:9,5]
q_ic3['AIC'] =q.iloc[9:17,5]
q_ic5['AIC'] =q.iloc[17:25,5]
q_ic7['AIC'] =q.iloc[25:33,5]
q_ic['1 Feature']=q_ic1['AIC'].reset_index(drop=True)
q_ic['3 Features']=q_ic3['AIC'].reset_index(drop=True)
q_ic['5 Features']=q_ic5['AIC'].reset_index(drop=True)
q_ic['7 Features']=q_ic7['AIC'].reset_index(drop=True)
q_ic = q_ic.T
q_ic.columns = q_ic.iloc[0,:]
q_ic = q_ic.iloc[1:,:]
print(q_ic)
color =['blue','orange','green','red','plum','brown','pink','gray']
y=['LinearRegression','ElasticNet','kNN','SVM','DecisionTree','RandomForest','AdaBoost','NeuralNet']
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)
xlabels=['1 Feaure','3 Features','5 Features','7 Features']
ax.set_xticks(np.arange(0,len(xlabels),1))
ax.set_xticklabels(xlabels, fontname = 'Arial', fontsize=15)
ax.set_xlabel('ML Model: Number of Features', fontname = 'Arial', fontsize=20)
ax.set_ylabel('AIC Score', fontname = 'Arial', fontsize=20)
plt.yticks(fontname = 'Arial', fontsize=15)


q_ic.plot(y=y,ax=ax,marker='.',markersize=12,color=color, lw=4)
#ax.set_title('TCR-pMHC Catch H-Bonds', fontname = 'Arial', fontsize=20)
#ax.set_xticklabels(['1 Feature','3 Features','5 Features','7 Features'], rotation=0)
ax.set_title('AIC for Quaternary Physiochemical Features', fontname = 'Arial', fontsize=20)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures_aic.png', bbox_inches='tight', dpi=300)

q_ic = pd.DataFrame()
q_ic1 = pd.DataFrame()
q_ic3 = pd.DataFrame()
q_ic5= pd.DataFrame()
q_ic7 = pd.DataFrame()
q_ic['Model']=q.iloc[1:9,0].reset_index(drop=True)
q_ic1['BIC'] = q.iloc[1:9,6]
q_ic3['BIC'] =q.iloc[9:17,6]
q_ic5['BIC'] =q.iloc[17:25,6]
q_ic7['BIC'] =q.iloc[25:33,6]
q_ic['1 Feature']=q_ic1['BIC'].reset_index(drop=True)
q_ic['3 Features']=q_ic3['BIC'].reset_index(drop=True)
q_ic['5 Features']=q_ic5['BIC'].reset_index(drop=True)
q_ic['7 Features']=q_ic7['BIC'].reset_index(drop=True)
q_ic = q_ic.T
q_ic.columns = q_ic.iloc[0,:]
q_ic = q_ic.iloc[1:,:]
print(q_ic)
y=['LinearRegression','ElasticNet','kNN','SVM','DecisionTree','RandomForest','AdaBoost','NeuralNet']
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)

ax.set_xlabel('ML Model: Number of Features', fontname = 'Arial', fontsize=20)
ax.set_ylabel('BIC Score', fontname = 'Arial', fontsize=20)
ax.set_xticks(np.arange(0,len(xlabels),1))
ax.set_xticklabels(xlabels, fontname = 'Arial', fontsize=15)
plt.yticks(fontname = 'Arial', fontsize=15)

q_ic.plot(y=y,ax=ax,marker='.',markersize=12,color=color, lw=4)
#ax.set_title('TCR-pMHC Catch H-Bonds', fontname = 'Arial', fontsize=20)
#ax.set_xticklabels(['1 Feature','3 Features','5 Features','7 Features'], rotation=0)
ax.set_title('BIC for Quaternary Physiochemical Features', fontname = 'Arial', fontsize=20)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures_bic.png', bbox_inches='tight', dpi=300)

###################################

#Primary Features#

###################################
q = pd.read_excel(r'/Users/zrollins/Documents/Documents/TCR_ML/InformationCriteria.xlsx', sheet_name='pfeatures', index_col = 0, header=None)


q_ic = pd.DataFrame()
q_ic1 = pd.DataFrame()
q_ic3 = pd.DataFrame()
q_ic5= pd.DataFrame()
q_ic7 = pd.DataFrame()
q_ic['Model']=q.iloc[1:9,0].reset_index(drop=True)
q_ic1['AIC'] = q.iloc[1:9,5]
q_ic3['AIC'] =q.iloc[9:17,5]
q_ic5['AIC'] =q.iloc[17:25,5]
q_ic7['AIC'] =q.iloc[25:33,5]
q_ic['1 Feature']=q_ic1['AIC'].reset_index(drop=True)
q_ic['3 Features']=q_ic3['AIC'].reset_index(drop=True)
q_ic['5 Features']=q_ic5['AIC'].reset_index(drop=True)
q_ic['7 Features']=q_ic7['AIC'].reset_index(drop=True)
q_ic = q_ic.T
q_ic.columns = q_ic.iloc[0,:]
q_ic = q_ic.iloc[1:,:]
print(q_ic)
y=['LinearRegression','ElasticNet','kNN','SVM','DecisionTree','RandomForest','AdaBoost','NeuralNet']
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)

ax.set_xlabel('ML Model: Number of Features', fontname = 'Arial', fontsize=20)
ax.set_ylabel('AIC Score', fontname = 'Arial', fontsize=20)
ax.set_xticks(np.arange(0,len(xlabels),1))
ax.set_xticklabels(xlabels, fontname = 'Arial', fontsize=15)
plt.yticks(fontname = 'Arial', fontsize=15)

q_ic.plot(y=y,ax=ax,marker='.',markersize=12,color=color, lw=4)
#ax.set_title('TCR-pMHC Catch H-Bonds', fontname = 'Arial', fontsize=20)
#ax.set_xticklabels(['1 Feature','3 Features','5 Features','7 Features'], rotation=0)
ax.set_title('AIC for Primary Physiochemical Features', fontname = 'Arial', fontsize=20)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures_aic.png', bbox_inches='tight', dpi=300)

q_ic = pd.DataFrame()
q_ic1 = pd.DataFrame()
q_ic3 = pd.DataFrame()
q_ic5= pd.DataFrame()
q_ic7 = pd.DataFrame()
q_ic['Model']=q.iloc[1:9,0].reset_index(drop=True)
q_ic1['BIC'] = q.iloc[1:9,6]
q_ic3['BIC'] =q.iloc[9:17,6]
q_ic5['BIC'] =q.iloc[17:25,6]
q_ic7['BIC'] =q.iloc[25:33,6]
q_ic['1 Feature']=q_ic1['BIC'].reset_index(drop=True)
q_ic['3 Features']=q_ic3['BIC'].reset_index(drop=True)
q_ic['5 Features']=q_ic5['BIC'].reset_index(drop=True)
q_ic['7 Features']=q_ic7['BIC'].reset_index(drop=True)
q_ic = q_ic.T
q_ic.columns = q_ic.iloc[0,:]
q_ic = q_ic.iloc[1:,:]
print(q_ic)
y=['LinearRegression','ElasticNet','kNN','SVM','DecisionTree','RandomForest','AdaBoost','NeuralNet']
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)

ax.set_xlabel('ML Model: Number of Features', fontname = 'Arial', fontsize=20)
ax.set_ylabel('BIC Score', fontname = 'Arial', fontsize=20)
ax.set_xticks(np.arange(0,len(xlabels),1))
ax.set_xticklabels(xlabels, fontname = 'Arial', fontsize=15)
plt.yticks(fontname = 'Arial', fontsize=15)

q_ic.plot(y=y,ax=ax,marker='.',markersize=12,color=color, lw=4)
#ax.set_title('TCR-pMHC Catch H-Bonds', fontname = 'Arial', fontsize=20)
#ax.set_xticklabels(['1 Feature','3 Features','5 Features','7 Features'], rotation=0)
ax.set_title('BIC for Primary Physiochemical Features', fontname = 'Arial', fontsize=20)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures_bic.png', bbox_inches='tight', dpi=300)



