#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 20:14:26 2021

@author: zrollins
"""

from sklearn import datasets
from sklearn.model_selection import learning_curve
import matplotlib.pyplot as plt
import numpy as np
from numpy import mean
from numpy import std
from numpy import absolute
import pandas as pd
import matplotlib.pyplot as plt
from numpy import arange
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import ElasticNet
from sklearn import linear_model
from sklearn import svm
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from mlxtend.feature_selection import ExhaustiveFeatureSelector as EFS
from collections import OrderedDict
import logging


data = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures/data_t_qfs.csv', index_col = 0, header=None)
data = data.rename(columns=data.iloc[0])
data = data.drop(data.index[0])



Xn = data.iloc[:,:-1]
X_names = list(Xn.columns)

###FEATURE SELECTION
enet = pd.DataFrame()
#enet['Total H-Bonds'] = data['Total H-Bonds']
enet['Total Contacts'] = data['Total Contacts']
#enet['Instant H-Bonds'] = data['Instant H-Bonds']
#enet['Instant Contacts'] = data['Instant Contacts']
#enet['RXN Distance'] = data['RXN Distance']
#enet['Max Frequency Distance'] = data['Max Frequency Distance']
#enet['SASA_TCR'] = data['SASA_TCR']
#enet['SASA_pMHC'] = data['SASA_pMHC']
#enet['RMSF_TCR'] = data['RMSF_TCR']
#enet['RMSF_pMHC'] = data['RMSF_pMHC']
#enet['GYR_TCR'] = data['GYR_TCR']
#enet['GYRx_TCR'] = data['GYRx_TCR']
#enet['GYRy_TCR'] = data['GYRy_TCR']
#enet['GYRz_TCR'] = data['GYRz_TCR']
#enet['GYR_pMHC'] = data['GYR_pMHC']
#enet['GYRx_pMHC'] = data['GYRx_pMHC']
#enet['GYRy_pMHC'] = data['GYRy_pMHC']
#enet['GYRz_pMHC'] = data['GYRz_pMHC']

#enet['CDR3 Distance'] = data['CDR3 Distance']
#enet['SASA_pep'] = data['SASA_pep']
#enet['SASA_CDR1a'] = data['SASA_CDR1a']
#enet['SASA_CDR2a'] = data['SASA_CDR2a']
#enet['SASA_CDR3a'] = data['SASA_CDR3a']
#enet['SASA_CDR1b'] = data['SASA_CDR1b']
#enet['SASA_CDR2b'] = data['SASA_CDR2b']
#enet['SASA_CDR3b'] = data['SASA_CDR3b']
#enet['SASA_MHCa'] = data['SASA_MHCa']
#enet['SASA_MHCb'] = data['SASA_MHCb']
#enet['RMSF_pep'] = data['RMSF_pep']
#enet['RMSF_CDR1a'] = data['RMSF_CDR1a']
#enet['RMSF_CDR2a'] = data['RMSF_CDR2a']
#enet['RMSF_CDR3a'] = data['RMSF_CDR3a']
#enet['RMSF_CDR1b'] = data['RMSF_CDR1b']
#enet['RMSF_CDR2b'] = data['RMSF_CDR2b']
#enet['RMSF_CDR3b'] = data['RMSF_CDR3b']
#enet['RMSF_MHCa'] = data['RMSF_MHCa']
#enet['RMSF_MHCb'] = data['RMSF_MHCb']
#enet['GYR_pep'] = data['GYR_pep']
#enet['GYRx_pep'] = data['GYRx_pep']
#enet['GYRy_pep'] = data['GYRy_pep']
#enet['GYRz_pep'] = data['GYRz_pep']
#enet['GYR_CDR1a'] = data['GYR_CDR1a']
#enet['GYRx_CDR1a'] = data['GYRx_CDR1a']
#enet['GYRy_CDR1a'] = data['GYRy_CDR1a']
#enet['GYRz_CDR1a'] = data['GYRz_CDR1a']
#enet['GYR_CDR2a'] = data['GYR_CDR2a']
#enet['GYRx_CDR2a'] = data['GYRx_CDR2a']
#enet['GYRy_CDR2a'] = data['GYRy_CDR2a']
#enet['GYRz_CDR2a'] = data['GYRz_CDR2a']
#enet['GYR_CDR3a'] = data['GYR_CDR3a']
#enet['GYRx_CDR3a'] = data['GYRx_CDR3a']
#enet['GYRy_CDR3a'] = data['GYRy_CDR3a']
#enet['GYRz_CDR3a'] = data['GYRz_CDR3a']
#enet['GYR_CDR1b'] = data['GYR_CDR1b']
#enet['GYRx_CDR1b'] = data['GYRx_CDR1b']
#enet['GYRy_CDR1b'] = data['GYRy_CDR1b']
#enet['GYRz_CDR1b'] = data['GYRz_CDR1b']
#enet['GYR_CDR2b'] = data['GYR_CDR2b']
#enet['GYRx_CDR2b'] = data['GYRx_CDR2b']
#enet['GYRy_CDR2b'] = data['GYRy_CDR2b']
#enet['GYRz_CDR2b'] = data['GYRz_CDR2b']
#enet['GYR_CDR3b'] = data['GYR_CDR3b']
#enet['GYRx_CDR3b'] = data['GYRx_CDR3b']
#enet['GYRy_CDR3b'] = data['GYRy_CDR3b']
#enet['GYRz_CDR3b'] = data['GYRz_CDR3b']
#enet['GYR_MHCa'] = data['GYR_MHCa']
#enet['GYRx_MHCa'] = data['GYRx_MHCa']
#enet['GYRy_MHCa'] = data['GYRy_MHCa']
#enet['GYRz_MHCa'] = data['GYRz_MHCa']
#enet['GYR_MHCb'] = data['GYR_MHCb']
#enet['GYRx_MHCb'] = data['GYRx_MHCb']
#enet['GYRy_MHCb'] = data['GYRy_MHCb']
#enet['GYRz_MHCb'] = data['GYRz_MHCb']
#enet['Instant Contacts CDR1a_MHCa']=data['Instant Contacts CDR1a_MHCa']
#enet['Instant Contacts CDR1a_pep']=data['Instant Contacts CDR1a_pep']
#enet['Total Contacts CDR1a_MHCa']=data['Total Contacts CDR1a_MHCa']
#enet['Total Contacts CDR1a_pep']=data['Total Contacts CDR1a_pep']
#enet['Instant Contacts CDR2a_MHCa']=data['Instant Contacts CDR2a_MHCa']
#enet['Instant Contacts CDR2a_pep']=data['Instant Contacts CDR2a_pep'
#enet['Total Contacts CDR2a_MHCa']=data['Total Contacts CDR2a_MHCa']
#enet['Total Contacts CDR2a_pep']=data['Total Contacts CDR2a_pep']
#enet['Instant Contacts CDR3a_pep']=data['Instant Contacts CDR3a_pep']
#enet['Instant Contacts CDR3a_MHCa']=data['Instant Contacts CDR3a_MHCa']
#enet['Total Contacts CDR3a_pep']=data['Total Contacts CDR3a_pep']
#enet['Total Contacts CDR3a_MHCa']=data['Total Contacts CDR3a_MHCa']
#enet['Instant Contacts CDR1b_MHCb']=data['Instant Contacts CDR1b_MHCb']
#enet['Instant Contacts CDR1b_pep']=data['Instant Contacts CDR1b_pep']
#enet['Total Contacts CDR1b_MHCb']=data['Total Contacts CDR1b_MHCb']
#enet['Total Contacts CDR1b_pep']=data['Total Contacts CDR1b_pep']
#enet['Instant Contacts CDR2b_MHCb']=data['Instant Contacts CDR2b_MHCb']
#enet['Instant Contacts CDR2b_pep']=data['Instant Contacts CDR2b_pep']
#enet['Total Contacts CDR2b_MHCb']=data['Total Contacts CDR2b_MHCb']
#enet['Total Contacts CDR2b_pep']=data['Total Contacts CDR2b_pep']
#enet['Instant Contacts CDR3b_pep']=data['Instant Contacts CDR3b_pep']
#enet['Instant Contacts CDR3b_MHCb']=data['Instant Contacts CDR3b_MHCb']
#enet['Total Contacts CDR3b_pep']=data['Total Contacts CDR3b_pep']
#enet['Total Contacts CDR3b_MHCb']=data['Total Contacts CDR3b_MHCb']
#enet['Instant H-Bonds CDR1a_MHCa']=data['Instant H-Bonds CDR1a_MHCa']
#enet['Instant H-Bonds CDR1a_pep']=data['Instant H-Bonds CDR1a_pep']
#enet['Total H-Bonds CDR1a_MHCa']=data['Total H-Bonds CDR1a_MHCa']
#enet['Total H-Bonds CDR1a_pep']=data['Total H-Bonds CDR1a_pep']
#enet['Instant H-Bonds CDR2a_MHCa']=data['Instant H-Bonds CDR2a_MHCa']
#enet['Instant H-Bonds CDR2a_pep']=data['Instant H-Bonds CDR2a_pep']
#enet['Total H-Bonds CDR2a_MHCa']=data['Total H-Bonds CDR2a_MHCa']
#enet['Total H-Bonds CDR2a_pep']=data['Total H-Bonds CDR2a_pep']
#enet['Instant H-Bonds CDR3a_pep']=data['Instant H-Bonds CDR3a_pep']
#enet['Instant H-Bonds CDR3a_MHCa']=data['Instant H-Bonds CDR3a_MHCa']
#enet['Total H-Bonds CDR3a_pep']=data['Total H-Bonds CDR3a_pep']
#enet['Total H-Bonds CDR3a_MHCa']=data['Total H-Bonds CDR3a_MHCa']
#enet['Instant H-Bonds CDR1b_MHCb']=data['Instant H-Bonds CDR1b_MHCb']
#enet['Instant H-Bonds CDR1b_pep']=data['Instant H-Bonds CDR1b_pep']
#enet['Total H-Bonds CDR1b_MHCb']=data['Total H-Bonds CDR1b_MHCb']
#enet['Total H-Bonds CDR1b_pep']=data['Total H-Bonds CDR1b_pep']
#enet['Instant H-Bonds CDR2b_MHCb']=data['Instant H-Bonds CDR2b_MHCb']
#enet['Instant H-Bonds CDR2b_pep']=data['Instant H-Bonds CDR2b_pep']
#enet['Total H-Bonds CDR2b_MHCb']=data['Total H-Bonds CDR2b_MHCb']
#enet['Total H-Bonds CDR2b_pep']=data['Total H-Bonds CDR2b_pep']
#enet['Instant H-Bonds CDR3b_pep']=data['Instant H-Bonds CDR3b_pep']
#enet['Instant H-Bonds CDR3b_MHCb']=data['Instant H-Bonds CDR3b_MHCb']
#enet['Total H-Bonds CDR3b_pep']=data['Total H-Bonds CDR3b_pep']
#enet['Total H-Bonds CDR3b_MHCb']=data['Total H-Bonds CDR3b_MHCb']

enet['Bond Lifetime'] = data['Bond Lifetime']

e = enet.values
#e = data.values
e = e[:,:].astype(float)

print(Xn)
X, y = e[:, :-1], e[:, -1]

cv = RepeatedKFold(n_splits=3, n_repeats=3, random_state=1)

svm = svm.SVR('linear',
              3,
              10,
              100)


# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validation:

train_sizes, train_scores, test_scores, fit_times, _ = learning_curve(svm, X, y,scoring='neg_median_absolute_error', return_times=True, cv=cv)
#train_sizes=np.arange(0, train_sizes, 1)

param_range = np.arange(0, len(train_scores), 1)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)
fit_times_mean = np.mean(fit_times, axis=1)
fit_times_std = np.std(fit_times, axis=1)


_, axes = plt.subplots(1, 2, figsize=(20, 5))
 # Plot learning curve
axes[0].fill_between(train_sizes, train_scores_mean - train_scores_std,
                         train_scores_mean + train_scores_std, alpha=0.1,
                         color="r")
axes[0].fill_between(train_sizes, test_scores_mean - test_scores_std,
                         test_scores_mean + test_scores_std, alpha=0.1,
                         color="g")
axes[0].plot(train_sizes, train_scores_mean, 'o-', color="r",
                 label="Training score")
axes[0].plot(train_sizes, test_scores_mean, 'o-', color="g",
                 label="Cross-validation score")
axes[0].legend(loc="best")
axes[0].set_xlabel("Training Examples")
axes[0].set_ylabel("Negative Mean Absolute Error")
axes[0].set_title("Learning Curve: (Quaternary Feature SVM)")
# Plot n_samples vs fit_times
axes[1].grid()
axes[1].plot(train_sizes, fit_times_mean, 'o-')
axes[1].fill_between(train_sizes, fit_times_mean - fit_times_std,
                         fit_times_mean + fit_times_std, alpha=0.1)
axes[1].set_xlabel("Training Examples")
axes[1].set_ylabel("Fit Time (s)")
axes[1].set_title("Scalability of the Model: (Quaternary SVM)")
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/qfeatures/learningcurve.png', dpi=300, bbox_inches='tight')