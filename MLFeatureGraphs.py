#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 12:33:15 2021

@author: zrollins
"""

import numpy as np
from numpy import mean
from numpy import std
from numpy import absolute
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from numpy import arange






pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 1000)

q1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/qfeatures/pred_t_1fs_best_model.csv', index_col = 0, header=None)
q2 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/qfeatures/pred_t_2fs_best_model.csv', index_col = 0, header=None)
q3 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/qfeatures/pred_t_3fs_best_model.csv', index_col = 0, header=None)
q5 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/qfeatures/pred_t_5fs_best_model.csv', index_col = 0, header=None)
q7 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/qfeatures/pred_t_7fs_best_model.csv', index_col = 0, header=None)

p1 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/pfeatures/pred_t_1fs_best_model.csv', index_col = 0, header=None)
p2 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/pfeatures/pred_t_2fs_best_model.csv', index_col = 0, header=None)
p3 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/pfeatures/pred_t_3fs_best_model.csv', index_col = 0, header=None)
p5 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/pfeatures/pred_t_5fs_best_model.csv', index_col = 0, header=None)
p7 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/pfeatures/pred_t_7fs_best_model.csv', index_col = 0, header=None)

print(q1)
q_bic = pd.DataFrame()
q_bic['Models'] = q1[1]
q_bic['1 Feature(s)'] =q1[2]
q_bic['2 Feature(s)'] =q2[2]
q_bic['3 Feature(s)'] =q3[2]
q_bic['5 Feature(s)'] =q5[2]
q_bic['7 Feature(s)'] =q7[2]

p_bic = pd.DataFrame()
p_bic['Models'] = p1[1]
p_bic['1 Feature(s)'] =p1[2]
p_bic['2 Feature(s)'] =p2[2]
p_bic['3 Feature(s)'] =p3[2]
p_bic['5 Feature(s)'] =p5[2]
p_bic['7 Feature(s)'] =p7[2]

q_p1 = pd.DataFrame()
q_p1['Models'] = q1[1]
q_p1['1 Feature(s) Mean'] =q1[2]
q_p1['2 Feature(s) Mean'] =q2[2]
q_p1['3 Feature(s) Mean'] =q3[2]
q_p1['5 Feature(s) Mean'] =q5[2]
q_p1['7 Feature(s) Mean'] =q7[2]
q_p2 = pd.DataFrame()
q_p2['Models'] = q1[1]
q_p2['1 Feature(s) std'] =q1[3]
q_p2['2 Feature(s) std'] =q2[3]
q_p2['3 Feature(s) std'] =q3[3]
q_p2['5 Feature(s) std'] =q5[3]
q_p2['7 Feature(s) std'] =q7[3]

p_p1 = pd.DataFrame()
p_p1['Models'] = p1[1]
p_p1['1 Feature(s) Mean'] =p1[2]
p_p1['2 Feature(s) Mean'] =p2[2]
p_p1['3 Feature(s) Mean'] =p3[2]
p_p1['5 Feature(s) Mean'] =p5[2]
p_p1['7 Feature(s) Mean'] =p7[2]
p_p2 = pd.DataFrame()
p_p2['Models'] = p1[1]
p_p2['1 Feature(s) std'] =p1[3]
p_p2['2 Feature(s) std'] =p2[3]
p_p2['3 Feature(s) std'] =p3[3]
p_p2['5 Feature(s) std'] =p5[3]
p_p2['7 Feature(s) std'] =p7[3]


#plot ML Model Performance
q_p1 = q_p1.iloc[1:,:].T
q_p2 = q_p2.iloc[1:,:].T
q_p1.columns =q_p1.iloc[0]
q_p2.columns =q_p2.iloc[0]
q_p1 = q_p1.iloc[1:,:].astype(float).div(5444).round(10).subtract(1).multiply(-100)
q_p2 = q_p2.iloc[1:,:].astype(float).div(5444).div(3).round(10).multiply(-100)

print(q_p1)
q_1f = pd.DataFrame()
q_1f_std =pd.DataFrame()
q_1f = q_p1.iloc[0,:]
q_1f_std = q_p2.iloc[0,:]
q_2f = pd.DataFrame()
q_2f_std =pd.DataFrame()
q_2f = q_p1.iloc[1,:]
q_2f_std = q_p2.iloc[1,:]
q_3f = pd.DataFrame()
q_3f_std =pd.DataFrame()
q_3f = q_p1.iloc[2,:]
q_3f_std = q_p2.iloc[2,:]
q_5f = pd.DataFrame()
q_5f_std =pd.DataFrame()
q_5f = q_p1.iloc[3,:]
q_5f_std = q_p2.iloc[3,:]
q_7f = pd.DataFrame()
q_7f_std =pd.DataFrame()
q_7f = q_p1.iloc[4,:]
q_7f_std = q_p2.iloc[4,:]
print(q_7f)
labels = ['Linear Regression','Elastic Net','kNN','SVM','Decision Tree','Random Forest','AdaBoost','Neural Net']
color =['blue','orange','green','red','plum','brown','pink','gray']
x = np.arange(len(labels))
width=0.2
fig, ax = plt.subplots(figsize=(20,10))
rects2 = ax.bar(x-1.5*width,q_1f,yerr=q_1f_std,hatch='*',width=width,color=color,label='Glycosylated')
rects2 = ax.bar(x-0.75*width,q_2f,yerr=q_2f_std,hatch='xx',width=width,color=color,label='Glycosylated')
rects3 = ax.bar(x-width/2,q_3f,yerr=q_3f_std,hatch='o',width=width,color=color,label='Deglycosylated')
rects3 = ax.bar(x+width/2,q_5f,yerr=q_5f_std,hatch='-',width=width,color=color,label='Deglycosylated')
rects3 = ax.bar(x+1.5*width,q_7f,yerr=q_7f_std,hatch='x',width=width,color=color,label='Deglycosylated')
leg_qf1=mpatches.Patch(hatch='*',label='1 Feature',facecolor='white')
leg_qf2=mpatches.Patch(hatch='xx',label='2 Features',facecolor='white')
leg_qf3=mpatches.Patch(hatch='o',label='3 Features',facecolor='white')
leg_qf5=mpatches.Patch(hatch='-',label='5 Features',facecolor='white')
leg_qf7=mpatches.Patch(hatch='x',label='7 Features',facecolor='white')
plt.legend(handles=[leg_qf1,leg_qf2,leg_qf3,leg_qf5,leg_qf7],fontsize='xx-large')
plt.xticks(range(0,len(labels)), labels=labels, fontsize=15)
plt.xlabel('ML Models: Number of Features',fontsize=20)
plt.ylabel('Mean Accuracy (%)',fontsize=20)
plt.title('Quaternary Feature Bond Lifetime Prediction',fontsize=25)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/qfeatures.png', dpi=300, bbox_inches='tight')
plt.show()


p_p1 = p_p1.iloc[1:,:].T
p_p2 = p_p2.iloc[1:,:].T
p_p1.columns =p_p1.iloc[0]
p_p2.columns =p_p2.iloc[0]
p_p1 = p_p1.iloc[1:,:].astype(float).div(5444).round(10).subtract(1).multiply(-100)
p_p2 = p_p2.iloc[1:,:].astype(float).div(5444).div(3).round(10).multiply(-100)


print(p_p1)
p_1f = pd.DataFrame()
p_1f_std =pd.DataFrame()
p_1f = p_p1.iloc[0,:]
p_1f_std = p_p2.iloc[0,:]
p_2f = pd.DataFrame()
p_2f_std =pd.DataFrame()
p_2f = p_p1.iloc[1,:]
p_2f_std = p_p2.iloc[1,:]
p_3f = pd.DataFrame()
p_3f_std =pd.DataFrame()
p_3f = p_p1.iloc[2,:]
p_3f_std = p_p2.iloc[2,:]
p_5f = pd.DataFrame()
p_5f_std =pd.DataFrame()
p_5f = p_p1.iloc[3,:]
p_5f_std = p_p2.iloc[3,:]
p_7f = pd.DataFrame()
p_7f_std =pd.DataFrame()
p_7f = p_p1.iloc[4,:]
QZA;hrdzk3Q6[pUc
             p_7f_std = p_p2.iloc[4,:]
print(p_7f)
labels = ['Linear Regression','Elastic Net','kNN','SVM','Decision Tree','Random Forest','AdaBoost','Neural Net']
x = np.arange(len(labels))
width=0.2
fig, ax = plt.subplots(figsize=(20,10))
rects2 = ax.bar(x-1.5*width,p_1f,yerr=p_1f_std,hatch='*',width=width,color=color,label='Glycosylated')
rects2 = ax.bar(x-0.75*width,p_2f,yerr=p_2f_std,hatch='xx',width=width,color=color,label='Glycosylated')
rects3 = ax.bar(x-width/2,p_3f,yerr=p_3f_std,hatch='o',width=width,color=color,label='Deglycosylated')
rects3 = ax.bar(x+width/2,p_5f,yerr=p_5f_std,hatch='-',width=width,color=color,label='Deglycosylated')
rects3 = ax.bar(x+1.5*width,p_7f,yerr=p_7f_std,hatch='x',width=width,color=color,label='Deglycosylated')
leg_qf1=mpatches.Patch(hatch='*',label='1 Feature',facecolor='white')
leg_qf2=mpatches.Patch(hatch='xx',label='2 Feature',facecolor='white')
leg_qf3=mpatches.Patch(hatch='o',label='3 Features',facecolor='white')
leg_qf5=mpatches.Patch(hatch='-',label='5 Features',facecolor='white')
leg_qf7=mpatches.Patch(hatch='x',label='7 Features',facecolor='white')
plt.legend(handles=[leg_qf1,leg_qf2,leg_qf3,leg_qf5,leg_qf7],fontsize='xx-large')
plt.xticks(range(0,len(labels)), labels=labels, fontsize=15)
plt.xlabel('ML Models: Number of Features',fontsize=20)
plt.ylabel('Mean Accuracy (%)',fontsize=20)
plt.title('Primary Feature Bond Lifetime Prediction',fontsize=25)
plt.savefig('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/AI/median/pfeatures.png', dpi=300, bbox_inches='tight')
plt.show()

    