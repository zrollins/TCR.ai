#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 12:57:57 2021

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

q = pd.read_excel(r'/Users/zrollins/Documents/Documents/TCR_ML/InformationCriteria.xlsx', sheet_name='10_qfeatures',index_col = 0, header=None)
q_plt = pd.DataFrame()
q_plt['AVG'] = q[10]
q_plt['STD'] = q[11]
q_plt = q_plt.iloc[1:11,:]
q_plt = q_plt.iloc[::-1]
yerr = q_plt['STD'].values
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)

#colormap
data_color = [x / max(q_plt['AVG']) for x in q_plt['AVG']]
my_cmap = plt.cm.get_cmap('Greys')
colors = my_cmap(data_color)

q_plt.plot.barh(y='AVG', xerr=yerr, width=0.8,color=colors, rot=0,ax=ax)
ax.set_ylabel('Quaternary Feature', fontname='Arial', fontsize=15)
ax.set_xlabel('Mean Absolute Error (MAE)', fontname='Arial', fontsize=15)
ax.set_title('TCR-pMHC Top Ten Quaternary  Features', fontname = 'Arial', fontsize=20)
plt.xticks(fontsize=20, fontname='Arial')
#plt.yticks(fontsize=20,fontname='Arial')
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures_top10.png', bbox_inches='tight', dpi=300)


q = pd.read_excel(r'/Users/zrollins/Documents/Documents/TCR_ML/InformationCriteria.xlsx', sheet_name='10_pfeatures',index_col = 0, header=None)
q_plt = pd.DataFrame()
q_plt['AVG'] = q[10]
q_plt['STD'] = q[11]
q_plt = q_plt.iloc[1:11,:]
q_plt = q_plt.iloc[::-1]
yerr = q_plt['STD'].values
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)

#colormap
data_color = [x / max(q_plt['AVG']) for x in q_plt['AVG']]
my_cmap = plt.cm.get_cmap('Greys')
colors = my_cmap(data_color)

q_plt.plot.barh(y='AVG', xerr=yerr, width=0.8,color=colors, rot=0,ax=ax)
ax.set_ylabel('Primary Feature', fontname='Arial', fontsize=15)
ax.set_xlabel('Mean Absolute Error (MAE)', fontname='Arial', fontsize=15)
ax.set_title('TCR-pMHC Top Ten Primary  Features', fontname = 'Arial', fontsize=20)
plt.xticks(fontsize=20, fontname='Arial')
#plt.yticks(fontsize=20, fontname='Arial')
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures_top10.png', bbox_inches='tight', dpi=300)


print(q_plt)