#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 19:08:55 2021

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

q1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures/data_t_qfs.csv', index_col = 0, header=0)
p1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/instantless/data_t_pfs.csv', index_col = 0, header=0)

print(q1)
q = pd.DataFrame()
q['Bond Lifetime'] = q1['Bond Lifetime']
q['Total Contacts'] = q1['Total Contacts']
q['Total H-Bonds'] = q1['Total H-Bonds']
q['Instant Contacts'] = q1['Instant Contacts']
q['Instant H-Bonds'] = q1['Instant H-Bonds']
q['RXN Distance'] = q1['RXN Distance']
q['GYRx_pMHC'] = q1['GYRx_pMHC']
q['GYRx_TCR'] = q1['GYRx_TCR']
q['RMSF_TCR'] = q1['RMSF_TCR']
q['GYRy_pMHC'] = q1['GYRy_pMHC']
q['GYRz_TCR'] = q1['GYRz_TCR']
q = q.iloc[1:,:]

q.corr(method='pearson').to_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures/pearson_q.csv')
q.corr(method='spearman').to_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures/spearman_q.csv')
print(q.corr(method='pearson'))
print(q.corr(method='spearman'))

p = pd.DataFrame()
p['Bond Lifetime'] = p1['Bond Lifetime']
p['Total H-Bonds CDR2b_MHCb'] = p1['Total H-Bonds CDR2b_MHCb']
p['Total Contacts CDR2b_MHCb'] = p1['Total Contacts CDR2b_MHCb']
p['Total H-Bonds CDR3b_pep'] = p1['Total H-Bonds CDR3b_pep']
p['Total Contacts CDR3b_pep'] = p1['Total Contacts CDR3b_pep']
p['Total Contacts CDR1a_pep'] = p1['Total Contacts CDR1a_pep']
p['Total H-Bonds CDR1a_pep'] = p1['Total H-Bonds CDR1a_pep']
p['Total Contacts CDR3b_MHCb'] = p1['Total Contacts CDR3b_MHCb']
p['SASA_MHCa'] = p1['SASA_MHCa']
p['Total Contacts CDR3a_pep'] = p1['Total Contacts CDR3a_pep']
p['GYRz_MHCb'] = p1['GYRz_MHCb']
p = p.iloc[1:,:]

p.corr(method='pearson').to_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/pearson_p.csv')
p.corr(method='spearman').to_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/spearman_p.csv')
print(p.corr(method='pearson'))
print(p.corr(method='spearman'))