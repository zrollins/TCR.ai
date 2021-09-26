#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 18:14:55 2021

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

q = pd.read_excel(r'/Users/zrollins/Documents/Documents/TCR_ML/InformationCriteria.xlsx', sheet_name='time',index_col = 0, header=None)
q_plt = pd.DataFrame()
q_plt['tavg'] = q[7]
q_plt['tse'] = q[9]
q_plt = q_plt.iloc[1:18,:]
print(q_plt)
q_plt = q_plt.reindex(['MART1','L1','GVA','8S','6V','hCD9','HSV1gp3','ImrA','Mtub1','Mtub2','5D','6Y','6H','5H','4Y','3F','2P'])
q_plt = q_plt.iloc[::-1]
yerr = q_plt['tse'].values
fig, ax = plt.subplots(figsize=(10, 10))
plt.sca(ax)

#colormap
data_color = [x / max(q_plt['tavg']) for x in q_plt['tavg']]
my_cmap = plt.cm.get_cmap('cool')
colors = my_cmap(data_color)
color =['gray','gray','gray','gray','gray','gray','gray','black','black','black','black','black','black','black','black','black','black']

q_plt.plot.barh(y='tavg', xerr=yerr, width=0.8,color=color, rot=0,ax=ax)
ax.set_ylabel('Epitope', fontname='Arial', fontsize=20)
ax.set_xlabel('Bond Lifetime (ps)', fontname='Arial', fontsize=20)
ax.set_title('TCR-pMHC Bond Lifetime (500 pN)', fontname = 'Arial', fontsize=20)
ax.legend([])
plt.xticks(fontsize=15)
plt.yticks(fontsize=20)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/bond_lifetime.png', bbox_inches='tight', dpi=300)


