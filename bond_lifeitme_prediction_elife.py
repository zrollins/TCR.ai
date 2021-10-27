#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 07:10:49 2021

@author: zrollins
"""

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

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
logger.addHandler(logging.FileHandler('log.log', 'a'))
print = logger.info

data = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/instantless/data_t_pfs.csv', index_col = 0, header=None)
data = data.rename(columns=data.iloc[0])
data = data.drop(data.index[0])

Xn = data.iloc[:,:-1]
X_names = list(Xn.columns)

###FEATURE SELECTION
enet = pd.DataFrame()
#enet['Total H-Bonds'] = data['Total H-Bonds']
#enet['Total Contacts'] = data['Total Contacts']
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
enet['SASA_MHCa'] = data['SASA_MHCa']
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
enet['GYRz_MHCb'] = data['GYRz_MHCb']
#enet['Instant Contacts CDR1a_MHCa']=data['Instant Contacts CDR1a_MHCa']
#enet['Instant Contacts CDR1a_pep']=data['Instant Contacts CDR1a_pep']
#enet['Total Contacts CDR1a_MHCa']=data['Total Contacts CDR1a_MHCa']
enet['Total Contacts CDR1a_pep']=data['Total Contacts CDR1a_pep']
#enet['Instant Contacts CDR2a_MHCa']=data['Instant Contacts CDR2a_MHCa']
#enet['Instant Contacts CDR2a_pep']=data['Instant Contacts CDR2a_pep'
#enet['Total Contacts CDR2a_MHCa']=data['Total Contacts CDR2a_MHCa']
#enet['Total Contacts CDR2a_pep']=data['Total Contacts CDR2a_pep']
#enet['Instant Contacts CDR3a_pep']=data['Instant Contacts CDR3a_pep']
#enet['Instant Contacts CDR3a_MHCa']=data['Instant Contacts CDR3a_MHCa']
enet['Total Contacts CDR3a_pep']=data['Total Contacts CDR3a_pep']
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
enet['Total Contacts CDR3b_pep']=data['Total Contacts CDR3b_pep']
enet['Total Contacts CDR3b_MHCb']=data['Total Contacts CDR3b_MHCb']
#enet['Instant H-Bonds CDR1a_MHCa']=data['Instant H-Bonds CDR1a_MHCa']
#enet['Instant H-Bonds CDR1a_pep']=data['Instant H-Bonds CDR1a_pep']
#enet['Total H-Bonds CDR1a_MHCa']=data['Total H-Bonds CDR1a_MHCa']
enet['Total H-Bonds CDR1a_pep']=data['Total H-Bonds CDR1a_pep']
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
enet['Total H-Bonds CDR2b_MHCb']=data['Total H-Bonds CDR2b_MHCb']
#enet['Total H-Bonds CDR2b_pep']=data['Total H-Bonds CDR2b_pep']
#enet['Instant H-Bonds CDR3b_pep']=data['Instant H-Bonds CDR3b_pep']
#enet['Instant H-Bonds CDR3b_MHCb']=data['Instant H-Bonds CDR3b_MHCb']
enet['Total H-Bonds CDR3b_pep']=data['Total H-Bonds CDR3b_pep']
#enet['Total H-Bonds CDR3b_MHCb']=data['Total H-Bonds CDR3b_MHCb']

enet['Bond Lifetime'] = data['Bond Lifetime']

e = enet.values
#e = data.values
e = e[:,:].astype(float)
                         

print(Xn)
X, y = e[:, :-1], e[:, -1]

cv = RepeatedKFold(n_splits=3, n_repeats=3, random_state=1)

# define models
model_params = {
    'LinearRegression':{
        'model': linear_model.TweedieRegressor(),
        'params':{
            'power':[0,1,2,3],
            'alpha':[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
        }
    },
    'ElasticNet':{
        'model':ElasticNet(),
        'params':{
           'alpha':[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0],
           'l1_ratio':arange(0, 1, 0.01)

        }
    },
    'kNN':{
        'model':KNeighborsRegressor(),
        'params':{
            'n_neighbors':arange(1,10,1),
            'leaf_size':arange(1,30,1),
            'p':[1,2]
        }
    },
    'SVM':{
        'model': svm.SVR(),
        'params':{
            'C':[0.1, 1, 10, 100, 1000],
            'gamma': [10, 1, 0.1, 0.01, 0.001, 0.0001],
            'kernel':['rbf','linear'],
            'degree':[3,4,5]
        }
    },
    'DecisionTree':{
        'model':DecisionTreeRegressor(),
        'params':{
            'splitter':['best','random'],
            'max_depth' :arange(1,10,1) ,
            'max_features':['auto', 'sqrt','log2'],
            'max_leaf_nodes':arange(2,10,1)
        }
    },
    'RandomForest':{
        'model':RandomForestRegressor(),
        'params':{
            'n_estimators':arange(10,1000,50),
            'max_features':['auto', 'sqrt','log2']
        }
    },
    'AdaBoost':{
        'model':AdaBoostRegressor(),
        'params':{
            'n_estimators':arange(10,1000,50),
            'loss':['linear','square','exponential']
        }
   },
    'NeuralNet':{
        'model':MLPRegressor(),
        'params':{
            'hidden_layer_sizes':arange(1,10,1),
            'activation':['relu','tanh','logistic','identity'],
            'alpha':[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
        }
    }
}

scores = []
cv_results = []


for model_name, mp in model_params.items():
    search =  GridSearchCV(mp['model'], mp['params'], scoring='neg_median_absolute_error', cv=cv,refit=False, n_jobs=-1)
    results = search.fit(X, y)
    cv_results.append({
        'model': model_name,
        'means': results.cv_results_['mean_test_score'],
        'stds': results.cv_results_['std_test_score'],
        'params': results.cv_results_['params'],
    })
    scores.append({
        'model': model_name,
        'best_score': results.best_score_,
        'best_stds':results.cv_results_['std_test_score'][results.best_index_],
        'best_params': results.best_params_,
    })
    
    
    
df = pd.DataFrame(scores,columns=['model','best_score','best_stds','best_params'])
df['best_score']=df['best_score']*-1
cv_results2 = pd.DataFrame(cv_results,columns=['model','means','stds','params'])
#colormap
data_color = [x / max(df['best_score']) for x in df['best_score']]
my_cmap = plt.cm.get_cmap('cool')
colors = my_cmap(data_color)
#plot
ax = df.plot(x='model',y='best_score', xerr='best_stds',kind='barh', color=colors)
ax.set_xlabel('Repeated 3-Fold Mean Test Score (Median Absolute Error)')
ax.set_ylabel('Machine Learning Model')
ax.set_title('TCR-pMHC Bond Lifetime Prediction (9 Primary Features)')
ax.get_legend().remove()
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/instantless/t_9fs.png', dpi=300, bbox_inches='tight')
#save data
df.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/instantless/t_9fs_best_model.csv')
cv_results2.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/pfeatures/instantless/t_9fs_cv.csv')

