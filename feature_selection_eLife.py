import numpy as np
from numpy import mean
from numpy import std
from numpy import absolute
import pandas as pd
import matplotlib.pyplot as plt
from numpy import arange
from sklearn.model_selection import train_test_split
#from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from mlxtend.feature_selection import ExhaustiveFeatureSelector as EFS
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
logger.addHandler(logging.FileHandler('log.log', 'a'))
print = logger.info


data = pd.read_csv('data_t_qfs.csv', index_col = 0, header=None)
data = data.rename(columns=data.iloc[0])
data = data.drop(data.index[0])
enet=pd.DataFrame()

#enet = data.iloc[:,146:]
enet['Total H-Bonds'] = data['Total H-Bonds']
enet['Total Contacts'] = data['Total Contacts']
enet['Instant H-Bonds'] = data['Instant H-Bonds']
enet['Instant Contacts'] = data['Instant Contacts']
enet['RXN Distance'] = data['RXN Distance']
#enet['Max Frequency Distance'] = data['Max Frequency Distance']
#enet['SASA_TCR'] = data['SASA_TCR']
#enet['SASA_pMHC'] = data['SASA_pMHC']
enet['RMSF_TCR'] = data['RMSF_TCR']
#enet['RMSF_pMHC'] = data['RMSF_pMHC']
#enet['GYR_TCR'] = data['GYR_TCR']
enet['GYRx_TCR'] = data['GYRx_TCR']
#enet['GYRy_TCR'] = data['GYRy_TCR']
enet['GYRz_TCR'] = data['GYRz_TCR']
#enet['GYR_pMHC'] = data['GYR_pMHC']
enet['GYRx_pMHC'] = data['GYRx_pMHC']
enet['GYRy_pMHC'] = data['GYRy_pMHC']
#enet['GYRz_pMHC'] = data['GYRz_pMHC']
enet['Bond Lifetime'] = data['Bond Lifetime']

Xn = enet.iloc[:,:-1]
X_names = list(Xn.columns)

e = enet.values
#e = data.values
print(Xn)
X, y = e[:, :-1], e[:, -1]
print(data)

# define model
model = ElasticNet()
# define model evaluation method
cv = RepeatedKFold(n_splits=3,
                   n_repeats=3,
                   random_state=1)
# define grid
grid = dict()
grid['exhaustivefeatureselector__estimator__alpha'] = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
grid['exhaustivefeatureselector__estimator__l1_ratio'] = arange(0, 1, 0.01)
# define search
efs = EFS(estimator=model,
           min_features=8,
           max_features=8,
           scoring='neg_mean_absolute_error',
           print_progress=False,
           clone_estimator=False,
           cv=cv,
           n_jobs=-1)
pipe = make_pipeline(efs, model)
search = GridSearchCV(estimator=pipe,
                      param_grid=grid,
                      scoring='neg_median_absolute_error',
                      cv=cv,
                      n_jobs=-1)
# perform the search & exhaustive feature selection
results = search.fit(X, y)
#summarize search
print(results)
print('MAE: %.3f' % results.best_score_)
print('Config: %s' % results.best_params_)
fidx = list(results.best_estimator_.steps[0][1].best_idx_)
features = [X_names[i] for i in fidx]
print('Best features:',features)
#print('Coeffs: %s' % results.best_estimator_.steps[0][1].coef_)
#print('Inter: %s' % results.best_estimator_.steps[0][1].intercept_)
df = pd.DataFrame(results.cv_results_)
df2 = pd.DataFrame(results.best_estimator_.steps[0][1].subsets_)
df3 =pd.DataFrame(features)

df.to_csv('cv_q8fs.csv')
df2.to_csv('grid_search_q8fs.csv')
df3.to_csv('features_q8fs.csv')
                                   