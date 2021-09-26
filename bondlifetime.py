#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 14:51:56 2021

@author: zrollins
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#bndx_time
t, t1, t2, bndx, bndx1, bndx2 = [],[],[],[],[],[]


##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/100/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            bndx.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/95/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            bndx1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/MART1/90/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            bndx2.append(float(cols[1]))
            
time1 = len(t)
time2 = len(t1) 
time3 = len(t2)

data_labels = {'MART1'}
data = pd.DataFrame(data_labels)
data['time1'] = time1
data['time2'] = time2
data['time3'] = time3
u = data.mean(0, skipna=True)
data['mean']=u
print(data)

#data.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/MART1/time.csv')

