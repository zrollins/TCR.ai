#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:39:57 2021

@author: zrollins
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 10:01:16 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#Residues

resi = ['Y','S','D','R','G','S','Q','S','F',
        'Y','S','N','G','D','K',
        'A','V','N','F','G','G','G','K','L','I','F',
        'Q','D','M','R','H','N','A',
        'Y','S','N','T','A','G','T','T',
        'A','S','S','L','S','F','G','T','E','A','F','F',
        'M','A','A','Q','T','T','K','H','K','W','E','A',
        'A','H','V','A','E','Q','L','R','A','Y','L','E',
        'G','T','C','V','E','W','L','R','R','Y','L','E',
        'N','G','K','E','T','L',
        'P','W','I','E','Q','E','G','P','E','Y','W','D',
        'G','E','T','R','K','V','K','A','H','S','Q','T',
        'H','R','V','D','L','G','T','L','R','G','Y','Y',
        'A','A','G','I','G','I','L','S','V']

#gyr
t, t1, t2, gyr_CDR1a, gyr_CDR1a1, gyr_CDR1a2, gyrx_CDR1a, gyrx_CDR1a1, gyrx_CDR1a2, gyry_CDR1a, gyry_CDR1a1, gyry_CDR1a2, gyrz_CDR1a, gyrz_CDR1a1, gyrz_CDR1a2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_CDR1a.append(float(cols[1]))
            gyrx_CDR1a.append(float(cols[2]))
            gyry_CDR1a.append(float(cols[3]))
            gyrz_CDR1a.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_CDR1a1.append(float(cols[1]))
            gyrx_CDR1a1.append(float(cols[2]))
            gyry_CDR1a1.append(float(cols[3]))
            gyrz_CDR1a1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_CDR1a2.append(float(cols[1]))
            gyrx_CDR1a2.append(float(cols[2]))
            gyry_CDR1a2.append(float(cols[3]))
            gyrz_CDR1a2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_CDR1a (1ns)': gyr_CDR1a}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_CDR1a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_CDR1a (2ns)': gyr_CDR1a1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_CDR1a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_CDR1a (3ns)': gyr_CDR1a2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_CDR1a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_CDR1a (2ns)'] = s1['gyr_CDR1a (2ns)']
s2 = s2.to_frame()
gyr['gyr_CDR1a (3ns)'] = s2['gyr_CDR1a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_CDR1a (1ns)','gyr_CDR1a (2ns)','gyr_CDR1a (3ns)']
gyr['gyr_CDR1a'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_CDR1a SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_CDR1a = gyr['gyr_CDR1a'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_CDR1a (1ns)': gyrx_CDR1a}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_CDR1a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_CDR1a (2ns)': gyrx_CDR1a1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_CDR1a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_CDR1a (3ns)': gyrx_CDR1a2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_CDR1a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_CDR1a (2ns)'] = s1['gyrx_CDR1a (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_CDR1a (3ns)'] = s2['gyrx_CDR1a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_CDR1a (1ns)','gyrx_CDR1a (2ns)','gyrx_CDR1a (3ns)']
gyrx['gyrx_CDR1a'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_CDR1a SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_CDR1a = gyrx['gyrx_CDR1a'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_CDR1a (1ns)': gyry_CDR1a}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_CDR1a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_CDR1a (2ns)': gyry_CDR1a1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_CDR1a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_CDR1a (3ns)': gyry_CDR1a2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_CDR1a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_CDR1a (2ns)'] = s1['gyry_CDR1a (2ns)']
s2 = s2.to_frame()
gyry['gyry_CDR1a (3ns)'] = s2['gyry_CDR1a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_CDR1a (1ns)','gyry_CDR1a (2ns)','gyry_CDR1a (3ns)']
gyry['gyry_CDR1a'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_CDR1a SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_CDR1a = gyry['gyry_CDR1a'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_CDR1a (1ns)': gyrz_CDR1a}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_CDR1a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_CDR1a (2ns)': gyrz_CDR1a1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_CDR1a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_CDR1a (3ns)': gyrz_CDR1a2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_CDR1a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_CDR1a (2ns)'] = s1['gyrz_CDR1a (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_CDR1a (3ns)'] = s2['gyrz_CDR1a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_CDR1a (1ns)','gyrz_CDR1a (2ns)','gyrz_CDR1a (3ns)']
gyrz['gyrz_CDR1a'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_CDR1a SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_CDR1a = gyrz['gyrz_CDR1a'].idxmax().right
#SASA
t, t1, t2, pep, pep1, pep2 = [],[],[],[],[],[]

#gyr
t, t1, t2, gyr_CDR2a, gyr_CDR2a1, gyr_CDR2a2, gyrx_CDR2a, gyrx_CDR2a1, gyrx_CDR2a2, gyry_CDR2a, gyry_CDR2a1, gyry_CDR2a2, gyrz_CDR2a, gyrz_CDR2a1, gyrz_CDR2a2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_CDR2a.append(float(cols[1]))
            gyrx_CDR2a.append(float(cols[2]))
            gyry_CDR2a.append(float(cols[3]))
            gyrz_CDR2a.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_CDR2a1.append(float(cols[1]))
            gyrx_CDR2a1.append(float(cols[2]))
            gyry_CDR2a1.append(float(cols[3]))
            gyrz_CDR2a1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_CDR2a2.append(float(cols[1]))
            gyrx_CDR2a2.append(float(cols[2]))
            gyry_CDR2a2.append(float(cols[3]))
            gyrz_CDR2a2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_CDR2a (1ns)': gyr_CDR2a}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_CDR2a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_CDR2a (2ns)': gyr_CDR2a1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_CDR2a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_CDR2a (3ns)': gyr_CDR2a2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_CDR2a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_CDR2a (2ns)'] = s1['gyr_CDR2a (2ns)']
s2 = s2.to_frame()
gyr['gyr_CDR2a (3ns)'] = s2['gyr_CDR2a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_CDR2a (1ns)','gyr_CDR2a (2ns)','gyr_CDR2a (3ns)']
gyr['gyr_CDR2a'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_CDR2a SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_CDR2a = gyr['gyr_CDR2a'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_CDR2a (1ns)': gyrx_CDR2a}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_CDR2a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_CDR2a (2ns)': gyrx_CDR2a1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_CDR2a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_CDR2a (3ns)': gyrx_CDR2a2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_CDR2a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_CDR2a (2ns)'] = s1['gyrx_CDR2a (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_CDR2a (3ns)'] = s2['gyrx_CDR2a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_CDR2a (1ns)','gyrx_CDR2a (2ns)','gyrx_CDR2a (3ns)']
gyrx['gyrx_CDR2a'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_CDR2a SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_CDR2a = gyrx['gyrx_CDR2a'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_CDR2a (1ns)': gyry_CDR2a}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_CDR2a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_CDR2a (2ns)': gyry_CDR2a1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_CDR2a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_CDR2a (3ns)': gyry_CDR2a2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_CDR2a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_CDR2a (2ns)'] = s1['gyry_CDR2a (2ns)']
s2 = s2.to_frame()
gyry['gyry_CDR2a (3ns)'] = s2['gyry_CDR2a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_CDR2a (1ns)','gyry_CDR2a (2ns)','gyry_CDR2a (3ns)']
gyry['gyry_CDR2a'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_CDR2a SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_CDR2a = gyry['gyry_CDR2a'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_CDR2a (1ns)': gyrz_CDR2a}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_CDR2a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_CDR2a (2ns)': gyrz_CDR2a1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_CDR2a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_CDR2a (3ns)': gyrz_CDR2a2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_CDR2a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_CDR2a (2ns)'] = s1['gyrz_CDR2a (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_CDR2a (3ns)'] = s2['gyrz_CDR2a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_CDR2a (1ns)','gyrz_CDR2a (2ns)','gyrz_CDR2a (3ns)']
gyrz['gyrz_CDR2a'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_CDR2a SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_CDR2a = gyrz['gyrz_CDR2a'].idxmax().right

#gyr
t, t1, t2, gyr_CDR3a, gyr_CDR3a1, gyr_CDR3a2, gyrx_CDR3a, gyrx_CDR3a1, gyrx_CDR3a2, gyry_CDR3a, gyry_CDR3a1, gyry_CDR3a2, gyrz_CDR3a, gyrz_CDR3a1, gyrz_CDR3a2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_CDR3a.append(float(cols[1]))
            gyrx_CDR3a.append(float(cols[2]))
            gyry_CDR3a.append(float(cols[3]))
            gyrz_CDR3a.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_CDR3a1.append(float(cols[1]))
            gyrx_CDR3a1.append(float(cols[2]))
            gyry_CDR3a1.append(float(cols[3]))
            gyrz_CDR3a1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_CDR3a2.append(float(cols[1]))
            gyrx_CDR3a2.append(float(cols[2]))
            gyry_CDR3a2.append(float(cols[3]))
            gyrz_CDR3a2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_CDR3a (1ns)': gyr_CDR3a}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_CDR3a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_CDR3a (2ns)': gyr_CDR3a1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_CDR3a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_CDR3a (3ns)': gyr_CDR3a2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_CDR3a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_CDR3a (2ns)'] = s1['gyr_CDR3a (2ns)']
s2 = s2.to_frame()
gyr['gyr_CDR3a (3ns)'] = s2['gyr_CDR3a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_CDR3a (1ns)','gyr_CDR3a (2ns)','gyr_CDR3a (3ns)']
gyr['gyr_CDR3a'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_CDR3a SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_CDR3a = gyr['gyr_CDR3a'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_CDR3a (1ns)': gyrx_CDR3a}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_CDR3a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_CDR3a (2ns)': gyrx_CDR3a1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_CDR3a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_CDR3a (3ns)': gyrx_CDR3a2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_CDR3a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_CDR3a (2ns)'] = s1['gyrx_CDR3a (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_CDR3a (3ns)'] = s2['gyrx_CDR3a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_CDR3a (1ns)','gyrx_CDR3a (2ns)','gyrx_CDR3a (3ns)']
gyrx['gyrx_CDR3a'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_CDR3a SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_CDR3a = gyrx['gyrx_CDR3a'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_CDR3a (1ns)': gyry_CDR3a}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_CDR3a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_CDR3a (2ns)': gyry_CDR3a1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_CDR3a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_CDR3a (3ns)': gyry_CDR3a2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_CDR3a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_CDR3a (2ns)'] = s1['gyry_CDR3a (2ns)']
s2 = s2.to_frame()
gyry['gyry_CDR3a (3ns)'] = s2['gyry_CDR3a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_CDR3a (1ns)','gyry_CDR3a (2ns)','gyry_CDR3a (3ns)']
gyry['gyry_CDR3a'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_CDR3a SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_CDR3a = gyry['gyry_CDR3a'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_CDR3a (1ns)': gyrz_CDR3a}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_CDR3a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_CDR3a (2ns)': gyrz_CDR3a1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_CDR3a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_CDR3a (3ns)': gyrz_CDR3a2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_CDR3a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_CDR3a (2ns)'] = s1['gyrz_CDR3a (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_CDR3a (3ns)'] = s2['gyrz_CDR3a (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_CDR3a (1ns)','gyrz_CDR3a (2ns)','gyrz_CDR3a (3ns)']
gyrz['gyrz_CDR3a'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_CDR3a SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_CDR3a = gyrz['gyrz_CDR3a'].idxmax().right

print(mx_gyr_CDR3a)
print(mx_gyrx_CDR3a)
print(mx_gyry_CDR3a)
print(mx_gyrz_CDR3a)

#gyr
t, t1, t2, gyr_CDR1b, gyr_CDR1b1, gyr_CDR1b2, gyrx_CDR1b, gyrx_CDR1b1, gyrx_CDR1b2, gyry_CDR1b, gyry_CDR1b1, gyry_CDR1b2, gyrz_CDR1b, gyrz_CDR1b1, gyrz_CDR1b2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_CDR1b.append(float(cols[1]))
            gyrx_CDR1b.append(float(cols[2]))
            gyry_CDR1b.append(float(cols[3]))
            gyrz_CDR1b.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_CDR1b1.append(float(cols[1]))
            gyrx_CDR1b1.append(float(cols[2]))
            gyry_CDR1b1.append(float(cols[3]))
            gyrz_CDR1b1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_CDR1b2.append(float(cols[1]))
            gyrx_CDR1b2.append(float(cols[2]))
            gyry_CDR1b2.append(float(cols[3]))
            gyrz_CDR1b2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_CDR1b (1ns)': gyr_CDR1b}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_CDR1b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_CDR1b (2ns)': gyr_CDR1b1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_CDR1b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_CDR1b (3ns)': gyr_CDR1b2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_CDR1b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_CDR1b (2ns)'] = s1['gyr_CDR1b (2ns)']
s2 = s2.to_frame()
gyr['gyr_CDR1b (3ns)'] = s2['gyr_CDR1b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_CDR1b (1ns)','gyr_CDR1b (2ns)','gyr_CDR1b (3ns)']
gyr['gyr_CDR1b'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_CDR1b SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_CDR1b = gyr['gyr_CDR1b'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_CDR1b (1ns)': gyrx_CDR1b}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_CDR1b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_CDR1b (2ns)': gyrx_CDR1b1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_CDR1b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_CDR1b (3ns)': gyrx_CDR1b2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_CDR1b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_CDR1b (2ns)'] = s1['gyrx_CDR1b (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_CDR1b (3ns)'] = s2['gyrx_CDR1b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_CDR1b (1ns)','gyrx_CDR1b (2ns)','gyrx_CDR1b (3ns)']
gyrx['gyrx_CDR1b'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_CDR1b SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_CDR1b = gyrx['gyrx_CDR1b'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_CDR1b (1ns)': gyry_CDR1b}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_CDR1b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_CDR1b (2ns)': gyry_CDR1b1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_CDR1b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_CDR1b (3ns)': gyry_CDR1b2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_CDR1b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_CDR1b (2ns)'] = s1['gyry_CDR1b (2ns)']
s2 = s2.to_frame()
gyry['gyry_CDR1b (3ns)'] = s2['gyry_CDR1b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_CDR1b (1ns)','gyry_CDR1b (2ns)','gyry_CDR1b (3ns)']
gyry['gyry_CDR1b'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_CDR1b SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_CDR1b = gyry['gyry_CDR1b'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_CDR1b (1ns)': gyrz_CDR1b}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_CDR1b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_CDR1b (2ns)': gyrz_CDR1b1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_CDR1b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_CDR1b (3ns)': gyrz_CDR1b2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_CDR1b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_CDR1b (2ns)'] = s1['gyrz_CDR1b (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_CDR1b (3ns)'] = s2['gyrz_CDR1b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_CDR1b (1ns)','gyrz_CDR1b (2ns)','gyrz_CDR1b (3ns)']
gyrz['gyrz_CDR1b'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_CDR1b SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_CDR1b = gyrz['gyrz_CDR1b'].idxmax().right

#gyr
t, t1, t2, gyr_CDR2b, gyr_CDR2b1, gyr_CDR2b2, gyrx_CDR2b, gyrx_CDR2b1, gyrx_CDR2b2, gyry_CDR2b, gyry_CDR2b1, gyry_CDR2b2, gyrz_CDR2b, gyrz_CDR2b1, gyrz_CDR2b2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_CDR2b.append(float(cols[1]))
            gyrx_CDR2b.append(float(cols[2]))
            gyry_CDR2b.append(float(cols[3]))
            gyrz_CDR2b.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_CDR2b1.append(float(cols[1]))
            gyrx_CDR2b1.append(float(cols[2]))
            gyry_CDR2b1.append(float(cols[3]))
            gyrz_CDR2b1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_CDR2b2.append(float(cols[1]))
            gyrx_CDR2b2.append(float(cols[2]))
            gyry_CDR2b2.append(float(cols[3]))
            gyrz_CDR2b2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_CDR2b (1ns)': gyr_CDR2b}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_CDR2b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_CDR2b (2ns)': gyr_CDR2b1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_CDR2b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_CDR2b (3ns)': gyr_CDR2b2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_CDR2b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_CDR2b (2ns)'] = s1['gyr_CDR2b (2ns)']
s2 = s2.to_frame()
gyr['gyr_CDR2b (3ns)'] = s2['gyr_CDR2b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_CDR2b (1ns)','gyr_CDR2b (2ns)','gyr_CDR2b (3ns)']
gyr['gyr_CDR2b'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_CDR2b SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_CDR2b = gyr['gyr_CDR2b'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_CDR2b (1ns)': gyrx_CDR2b}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_CDR2b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_CDR2b (2ns)': gyrx_CDR2b1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_CDR2b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_CDR2b (3ns)': gyrx_CDR2b2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_CDR2b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_CDR2b (2ns)'] = s1['gyrx_CDR2b (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_CDR2b (3ns)'] = s2['gyrx_CDR2b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_CDR2b (1ns)','gyrx_CDR2b (2ns)','gyrx_CDR2b (3ns)']
gyrx['gyrx_CDR2b'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_CDR2b SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_CDR2b = gyrx['gyrx_CDR2b'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_CDR2b (1ns)': gyry_CDR2b}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_CDR2b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_CDR2b (2ns)': gyry_CDR2b1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_CDR2b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_CDR2b (3ns)': gyry_CDR2b2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_CDR2b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_CDR2b (2ns)'] = s1['gyry_CDR2b (2ns)']
s2 = s2.to_frame()
gyry['gyry_CDR2b (3ns)'] = s2['gyry_CDR2b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_CDR2b (1ns)','gyry_CDR2b (2ns)','gyry_CDR2b (3ns)']
gyry['gyry_CDR2b'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_CDR2b SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_CDR2b = gyry['gyry_CDR2b'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_CDR2b (1ns)': gyrz_CDR2b}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_CDR2b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_CDR2b (2ns)': gyrz_CDR2b1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_CDR2b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_CDR2b (3ns)': gyrz_CDR2b2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_CDR2b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_CDR2b (2ns)'] = s1['gyrz_CDR2b (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_CDR2b (3ns)'] = s2['gyrz_CDR2b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_CDR2b (1ns)','gyrz_CDR2b (2ns)','gyrz_CDR2b (3ns)']
gyrz['gyrz_CDR2b'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_CDR2b SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_CDR2b = gyrz['gyrz_CDR2b'].idxmax().right

#gyr
t, t1, t2, gyr_CDR3b, gyr_CDR3b1, gyr_CDR3b2, gyrx_CDR3b, gyrx_CDR3b1, gyrx_CDR3b2, gyry_CDR3b, gyry_CDR3b1, gyry_CDR3b2, gyrz_CDR3b, gyrz_CDR3b1, gyrz_CDR3b2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_CDR3b.append(float(cols[1]))
            gyrx_CDR3b.append(float(cols[2]))
            gyry_CDR3b.append(float(cols[3]))
            gyrz_CDR3b.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_CDR3b1.append(float(cols[1]))
            gyrx_CDR3b1.append(float(cols[2]))
            gyry_CDR3b1.append(float(cols[3]))
            gyrz_CDR3b1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_CDR3b2.append(float(cols[1]))
            gyrx_CDR3b2.append(float(cols[2]))
            gyry_CDR3b2.append(float(cols[3]))
            gyrz_CDR3b2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_CDR3b (1ns)': gyr_CDR3b}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_CDR3b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_CDR3b (2ns)': gyr_CDR3b1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_CDR3b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_CDR3b (3ns)': gyr_CDR3b2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_CDR3b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_CDR3b (2ns)'] = s1['gyr_CDR3b (2ns)']
s2 = s2.to_frame()
gyr['gyr_CDR3b (3ns)'] = s2['gyr_CDR3b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_CDR3b (1ns)','gyr_CDR3b (2ns)','gyr_CDR3b (3ns)']
gyr['gyr_CDR3b'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_CDR3b SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_CDR3b = gyr['gyr_CDR3b'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_CDR3b (1ns)': gyrx_CDR3b}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_CDR3b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_CDR3b (2ns)': gyrx_CDR3b1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_CDR3b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_CDR3b (3ns)': gyrx_CDR3b2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_CDR3b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_CDR3b (2ns)'] = s1['gyrx_CDR3b (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_CDR3b (3ns)'] = s2['gyrx_CDR3b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_CDR3b (1ns)','gyrx_CDR3b (2ns)','gyrx_CDR3b (3ns)']
gyrx['gyrx_CDR3b'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_CDR3b SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_CDR3b = gyrx['gyrx_CDR3b'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_CDR3b (1ns)': gyry_CDR3b}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_CDR3b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_CDR3b (2ns)': gyry_CDR3b1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_CDR3b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_CDR3b (3ns)': gyry_CDR3b2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_CDR3b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_CDR3b (2ns)'] = s1['gyry_CDR3b (2ns)']
s2 = s2.to_frame()
gyry['gyry_CDR3b (3ns)'] = s2['gyry_CDR3b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_CDR3b (1ns)','gyry_CDR3b (2ns)','gyry_CDR3b (3ns)']
gyry['gyry_CDR3b'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_CDR3b SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_CDR3b = gyry['gyry_CDR3b'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_CDR3b (1ns)': gyrz_CDR3b}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_CDR3b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_CDR3b (2ns)': gyrz_CDR3b1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_CDR3b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_CDR3b (3ns)': gyrz_CDR3b2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_CDR3b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_CDR3b (2ns)'] = s1['gyrz_CDR3b (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_CDR3b (3ns)'] = s2['gyrz_CDR3b (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_CDR3b (1ns)','gyrz_CDR3b (2ns)','gyrz_CDR3b (3ns)']
gyrz['gyrz_CDR3b'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_CDR3b SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_CDR3b = gyrz['gyrz_CDR3b'].idxmax().right

#gyr
t, t1, t2, gyr_MHCa, gyr_MHCa1, gyr_MHCa2, gyrx_MHCa, gyrx_MHCa1, gyrx_MHCa2, gyry_MHCa, gyry_MHCa1, gyry_MHCa2, gyrz_MHCa, gyrz_MHCa1, gyrz_MHCa2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_MHCa.append(float(cols[1]))
            gyrx_MHCa.append(float(cols[2]))
            gyry_MHCa.append(float(cols[3]))
            gyrz_MHCa.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_MHCa1.append(float(cols[1]))
            gyrx_MHCa1.append(float(cols[2]))
            gyry_MHCa1.append(float(cols[3]))
            gyrz_MHCa1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_MHCa2.append(float(cols[1]))
            gyrx_MHCa2.append(float(cols[2]))
            gyry_MHCa2.append(float(cols[3]))
            gyrz_MHCa2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_MHCa (1ns)': gyr_MHCa}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_MHCa (2ns)': gyr_MHCa1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_MHCa (3ns)': gyr_MHCa2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_MHCa (2ns)'] = s1['gyr_MHCa (2ns)']
s2 = s2.to_frame()
gyr['gyr_MHCa (3ns)'] = s2['gyr_MHCa (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_MHCa (1ns)','gyr_MHCa (2ns)','gyr_MHCa (3ns)']
gyr['gyr_MHCa'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_MHCa SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_MHCa = gyr['gyr_MHCa'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_MHCa (1ns)': gyrx_MHCa}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_MHCa (2ns)': gyrx_MHCa1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_MHCa (3ns)': gyrx_MHCa2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_MHCa (2ns)'] = s1['gyrx_MHCa (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_MHCa (3ns)'] = s2['gyrx_MHCa (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_MHCa (1ns)','gyrx_MHCa (2ns)','gyrx_MHCa (3ns)']
gyrx['gyrx_MHCa'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_MHCa SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_MHCa = gyrx['gyrx_MHCa'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_MHCa (1ns)': gyry_MHCa}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_MHCa (2ns)': gyry_MHCa1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_MHCa (3ns)': gyry_MHCa2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_MHCa (2ns)'] = s1['gyry_MHCa (2ns)']
s2 = s2.to_frame()
gyry['gyry_MHCa (3ns)'] = s2['gyry_MHCa (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_MHCa (1ns)','gyry_MHCa (2ns)','gyry_MHCa (3ns)']
gyry['gyry_MHCa'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_MHCa SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_MHCa = gyry['gyry_MHCa'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_MHCa (1ns)': gyrz_MHCa}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_MHCa (2ns)': gyrz_MHCa1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_MHCa (3ns)': gyrz_MHCa2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_MHCa (2ns)'] = s1['gyrz_MHCa (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_MHCa (3ns)'] = s2['gyrz_MHCa (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_MHCa (1ns)','gyrz_MHCa (2ns)','gyrz_MHCa (3ns)']
gyrz['gyrz_MHCa'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_MHCa SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_MHCa = gyrz['gyrz_MHCa'].idxmax().right

#gyr
t, t1, t2, gyr_MHCb, gyr_MHCb1, gyr_MHCb2, gyrx_MHCb, gyrx_MHCb1, gyrx_MHCb2, gyry_MHCb, gyry_MHCb1, gyry_MHCb2, gyrz_MHCb, gyrz_MHCb1, gyrz_MHCb2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_MHCb.append(float(cols[1]))
            gyrx_MHCb.append(float(cols[2]))
            gyry_MHCb.append(float(cols[3]))
            gyrz_MHCb.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_MHCb1.append(float(cols[1]))
            gyrx_MHCb1.append(float(cols[2]))
            gyry_MHCb1.append(float(cols[3]))
            gyrz_MHCb1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_MHCb2.append(float(cols[1]))
            gyrx_MHCb2.append(float(cols[2]))
            gyry_MHCb2.append(float(cols[3]))
            gyrz_MHCb2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_MHCb (1ns)': gyr_MHCb}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_MHCb (2ns)': gyr_MHCb1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_MHCb (3ns)': gyr_MHCb2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_MHCb (2ns)'] = s1['gyr_MHCb (2ns)']
s2 = s2.to_frame()
gyr['gyr_MHCb (3ns)'] = s2['gyr_MHCb (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_MHCb (1ns)','gyr_MHCb (2ns)','gyr_MHCb (3ns)']
gyr['gyr_MHCb'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_MHCb SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_MHCb = gyr['gyr_MHCb'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_MHCb (1ns)': gyrx_MHCb}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_MHCb (2ns)': gyrx_MHCb1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_MHCb (3ns)': gyrx_MHCb2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_MHCb (2ns)'] = s1['gyrx_MHCb (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_MHCb (3ns)'] = s2['gyrx_MHCb (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_MHCb (1ns)','gyrx_MHCb (2ns)','gyrx_MHCb (3ns)']
gyrx['gyrx_MHCb'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_MHCb SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_MHCb = gyrx['gyrx_MHCb'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_MHCb (1ns)': gyry_MHCb}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_MHCb (2ns)': gyry_MHCb1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_MHCb (3ns)': gyry_MHCb2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_MHCb (2ns)'] = s1['gyry_MHCb (2ns)']
s2 = s2.to_frame()
gyry['gyry_MHCb (3ns)'] = s2['gyry_MHCb (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_MHCb (1ns)','gyry_MHCb (2ns)','gyry_MHCb (3ns)']
gyry['gyry_MHCb'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_MHCb SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_MHCb = gyry['gyry_MHCb'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_MHCb (1ns)': gyrz_MHCb}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_MHCb (2ns)': gyrz_MHCb1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_MHCb (3ns)': gyrz_MHCb2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_MHCb (2ns)'] = s1['gyrz_MHCb (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_MHCb (3ns)'] = s2['gyrz_MHCb (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_MHCb (1ns)','gyrz_MHCb (2ns)','gyrz_MHCb (3ns)']
gyrz['gyrz_MHCb'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_MHCb SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_MHCb = gyrz['gyrz_MHCb'].idxmax().right

#gyr
t, t1, t2, gyr_pep, gyr_pep1, gyr_pep2, gyrx_pep, gyrx_pep1, gyrx_pep2, gyry_pep, gyry_pep1, gyry_pep2, gyrz_pep, gyrz_pep1, gyrz_pep2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_pep.append(float(cols[1]))
            gyrx_pep.append(float(cols[2]))
            gyry_pep.append(float(cols[3]))
            gyrz_pep.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_pep1.append(float(cols[1]))
            gyrx_pep1.append(float(cols[2]))
            gyry_pep1.append(float(cols[3]))
            gyrz_pep1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_pep2.append(float(cols[1]))
            gyrx_pep2.append(float(cols[2]))
            gyry_pep2.append(float(cols[3]))
            gyrz_pep2.append(float(cols[4]))

#Bin gyr  
cut_bins=np.linspace(0,2,200).tolist()
            
gyr={'gyr_pep (1ns)': gyr_pep}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_pep (2ns)': gyr_pep1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_pep (3ns)': gyr_pep2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_pep (2ns)'] = s1['gyr_pep (2ns)']
s2 = s2.to_frame()
gyr['gyr_pep (3ns)'] = s2['gyr_pep (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_pep (1ns)','gyr_pep (2ns)','gyr_pep (3ns)']
gyr['gyr_pep'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_pep SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_pep = gyr['gyr_pep'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrx={'gyrx_pep (1ns)': gyrx_pep}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_pep (2ns)': gyrx_pep1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_pep (3ns)': gyrx_pep2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_pep (2ns)'] = s1['gyrx_pep (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_pep (3ns)'] = s2['gyrx_pep (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_pep (1ns)','gyrx_pep (2ns)','gyrx_pep (3ns)']
gyrx['gyrx_pep'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_pep SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_pep = gyrx['gyrx_pep'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,2,200).tolist()
            
gyry={'gyry_pep (1ns)': gyry_pep}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_pep (2ns)': gyry_pep1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_pep (3ns)': gyry_pep2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_pep (2ns)'] = s1['gyry_pep (2ns)']
s2 = s2.to_frame()
gyry['gyry_pep (3ns)'] = s2['gyry_pep (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_pep (1ns)','gyry_pep (2ns)','gyry_pep (3ns)']
gyry['gyry_pep'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_pep SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_pep = gyry['gyry_pep'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,2,200).tolist()
            
gyrz={'gyrz_pep (1ns)': gyrz_pep}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_pep (2ns)': gyrz_pep1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_pep (3ns)': gyrz_pep2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_pep (2ns)'] = s1['gyrz_pep (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_pep (3ns)'] = s2['gyrz_pep (3ns)']

l =np.linspace(0,2,200).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_pep (1ns)','gyrz_pep (2ns)','gyrz_pep (3ns)']
gyrz['gyrz_pep'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_pep SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_pep = gyrz['gyrz_pep'].idxmax().right


##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            pep2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'pep (1ns)': pep}
df = pd.DataFrame(sasa)
s=pd.cut(df['pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'pep (2ns)': pep1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'pep (3ns)': pep2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['pep (2ns)'] = s1['pep (2ns)']
s2 = s2.to_frame()
SASA['pep (3ns)'] = s2['pep (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['pep (1ns)','pep (2ns)','pep (3ns)']
SASA['pep'] = SASA[M1].astype(float).mean(axis=1)
SASA['pep SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_pep = SASA['pep'].idxmax().right

print(mx_pep)

#SASA
t, t1, t2, CDR1a, CDR1a1, CDR1a2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR1a.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR1a1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR1a2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'CDR1a (1ns)': CDR1a}
df = pd.DataFrame(sasa)
s=pd.cut(df['CDR1a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'CDR1a (2ns)': CDR1a1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['CDR1a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'CDR1a (3ns)': CDR1a2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['CDR1a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['CDR1a (2ns)'] = s1['CDR1a (2ns)']
s2 = s2.to_frame()
SASA['CDR1a (3ns)'] = s2['CDR1a (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['CDR1a (1ns)','CDR1a (2ns)','CDR1a (3ns)']
SASA['CDR1a'] = SASA[M1].astype(float).mean(axis=1)
SASA['CDR1a SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_CDR1a = SASA['CDR1a'].idxmax().right

print(mx_CDR1a)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:42:27 2021

@author: zrollins
"""

#SASA
t, t1, t2, CDR2a, CDR2a1, CDR2a2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR2a.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR2a1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR2a2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'CDR2a (1ns)': CDR2a}
df = pd.DataFrame(sasa)
s=pd.cut(df['CDR2a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'CDR2a (2ns)': CDR2a1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['CDR2a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'CDR2a (3ns)': CDR2a2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['CDR2a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['CDR2a (2ns)'] = s1['CDR2a (2ns)']
s2 = s2.to_frame()
SASA['CDR2a (3ns)'] = s2['CDR2a (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['CDR2a (1ns)','CDR2a (2ns)','CDR2a (3ns)']
SASA['CDR2a'] = SASA[M1].astype(float).mean(axis=1)
SASA['CDR2a SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_CDR2a = SASA['CDR2a'].idxmax().right

print(mx_CDR2a)

#SASA
t, t1, t2, CDR3a, CDR3a1, CDR3a2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR3a.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR3a1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR3a2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'CDR3a (1ns)': CDR3a}
df = pd.DataFrame(sasa)
s=pd.cut(df['CDR3a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'CDR3a (2ns)': CDR3a1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['CDR3a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'CDR3a (3ns)': CDR3a2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['CDR3a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['CDR3a (2ns)'] = s1['CDR3a (2ns)']
s2 = s2.to_frame()
SASA['CDR3a (3ns)'] = s2['CDR3a (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['CDR3a (1ns)','CDR3a (2ns)','CDR3a (3ns)']
SASA['CDR3a'] = SASA[M1].astype(float).mean(axis=1)
SASA['CDR3a SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_CDR3a = SASA['CDR3a'].idxmax().right

print(mx_CDR3a)

#SASA
t, t1, t2, CDR1b, CDR1b1, CDR1b2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR1b.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR1b1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR1b2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'CDR1b (1ns)': CDR1b}
df = pd.DataFrame(sasa)
s=pd.cut(df['CDR1b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'CDR1b (2ns)': CDR1b1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['CDR1b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'CDR1b (3ns)': CDR1b2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['CDR1b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['CDR1b (2ns)'] = s1['CDR1b (2ns)']
s2 = s2.to_frame()
SASA['CDR1b (3ns)'] = s2['CDR1b (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['CDR1b (1ns)','CDR1b (2ns)','CDR1b (3ns)']
SASA['CDR1b'] = SASA[M1].astype(float).mean(axis=1)
SASA['CDR1b SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_CDR1b = SASA['CDR1b'].idxmax().right

print(mx_CDR1b)


#SASA
t, t1, t2, CDR2b, CDR2b1, CDR2b2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR2b.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR2b1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR2b2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'CDR2b (1ns)': CDR2b}
df = pd.DataFrame(sasa)
s=pd.cut(df['CDR2b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'CDR2b (2ns)': CDR2b1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['CDR2b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'CDR2b (3ns)': CDR2b2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['CDR2b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['CDR2b (2ns)'] = s1['CDR2b (2ns)']
s2 = s2.to_frame()
SASA['CDR2b (3ns)'] = s2['CDR2b (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['CDR2b (1ns)','CDR2b (2ns)','CDR2b (3ns)']
SASA['CDR2b'] = SASA[M1].astype(float).mean(axis=1)
SASA['CDR2b SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_CDR2b = SASA['CDR2b'].idxmax().right

print(mx_CDR2b)

#SASA
t, t1, t2, CDR3b, CDR3b1, CDR3b2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR3b.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR3b1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR3b2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'CDR3b (1ns)': CDR3b}
df = pd.DataFrame(sasa)
s=pd.cut(df['CDR3b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'CDR3b (2ns)': CDR3b1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['CDR3b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'CDR3b (3ns)': CDR3b2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['CDR3b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['CDR3b (2ns)'] = s1['CDR3b (2ns)']
s2 = s2.to_frame()
SASA['CDR3b (3ns)'] = s2['CDR3b (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['CDR3b (1ns)','CDR3b (2ns)','CDR3b (3ns)']
SASA['CDR3b'] = SASA[M1].astype(float).mean(axis=1)
SASA['CDR3b SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_CDR3b = SASA['CDR3b'].idxmax().right

print(mx_CDR3b)

#SASA
t, t1, t2, MHCa, MHCa1, MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            MHCa.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            MHCa1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            MHCa2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'MHCa (1ns)': MHCa}
df = pd.DataFrame(sasa)
s=pd.cut(df['MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'MHCa (2ns)': MHCa1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'MHCa (3ns)': MHCa2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['MHCa (2ns)'] = s1['MHCa (2ns)']
s2 = s2.to_frame()
SASA['MHCa (3ns)'] = s2['MHCa (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['MHCa (1ns)','MHCa (2ns)','MHCa (3ns)']
SASA['MHCa'] = SASA[M1].astype(float).mean(axis=1)
SASA['MHCa SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_MHCa = SASA['MHCa'].idxmax().right

print(mx_MHCa)

#SASA
t, t1, t2, MHCb, MHCb1, MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            MHCb.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            MHCb1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            MHCb2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(1,100,1000).tolist()
            
sasa={'MHCb (1ns)': MHCb}
df = pd.DataFrame(sasa)
s=pd.cut(df['MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'MHCb (2ns)': MHCb1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'MHCb (3ns)': MHCb2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['MHCb (2ns)'] = s1['MHCb (2ns)']
s2 = s2.to_frame()
SASA['MHCb (3ns)'] = s2['MHCb (3ns)']

l =np.linspace(1,100,1000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['MHCb (1ns)','MHCb (2ns)','MHCb (3ns)']
SASA['MHCb'] = SASA[M1].astype(float).mean(axis=1)
SASA['MHCb SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_MHCb = SASA['MHCb'].idxmax().right

print(mx_MHCb)



#rmsf
t, t1, t2, pep, pep1, pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            pep2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'pep (1ns)': pep}
df = pd.DataFrame(rmsf)
s=pd.cut(df['pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'pep (2ns)': pep1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'pep (3ns)': pep2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['pep (2ns)'] = s1['pep (2ns)']
s2 = s2.to_frame()
rmsf['pep (3ns)'] = s2['pep (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['pep (1ns)','pep (2ns)','pep (3ns)']
rmsf['pep'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['pep SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_pep = rmsf['pep'].idxmax().right

print(mx_rmsf_pep)

#rmsf
t, t1, t2, CDR1a, CDR1a1, CDR1a2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR1a.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR1a1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_CDR1a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR1a2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'CDR1a (1ns)': CDR1a}
df = pd.DataFrame(rmsf)
s=pd.cut(df['CDR1a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'CDR1a (2ns)': CDR1a1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['CDR1a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'CDR1a (3ns)': CDR1a2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['CDR1a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['CDR1a (2ns)'] = s1['CDR1a (2ns)']
s2 = s2.to_frame()
rmsf['CDR1a (3ns)'] = s2['CDR1a (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['CDR1a (1ns)','CDR1a (2ns)','CDR1a (3ns)']
rmsf['CDR1a'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['CDR1a SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_CDR1a = rmsf['CDR1a'].idxmax().right

print(mx_rmsf_CDR1a)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:42:27 2021

@author: zrollins
"""

#rmsf
t, t1, t2, CDR2a, CDR2a1, CDR2a2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR2a.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR2a1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_CDR2a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR2a2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'CDR2a (1ns)': CDR2a}
df = pd.DataFrame(rmsf)
s=pd.cut(df['CDR2a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'CDR2a (2ns)': CDR2a1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['CDR2a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'CDR2a (3ns)': CDR2a2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['CDR2a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['CDR2a (2ns)'] = s1['CDR2a (2ns)']
s2 = s2.to_frame()
rmsf['CDR2a (3ns)'] = s2['CDR2a (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['CDR2a (1ns)','CDR2a (2ns)','CDR2a (3ns)']
rmsf['CDR2a'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['CDR2a SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_CDR2a = rmsf['CDR2a'].idxmax().right

print(mx_rmsf_CDR2a)

#rmsf
t, t1, t2, CDR3a, CDR3a1, CDR3a2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR3a.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR3a1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_CDR3a.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR3a2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'CDR3a (1ns)': CDR3a}
df = pd.DataFrame(rmsf)
s=pd.cut(df['CDR3a (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'CDR3a (2ns)': CDR3a1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['CDR3a (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'CDR3a (3ns)': CDR3a2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['CDR3a (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['CDR3a (2ns)'] = s1['CDR3a (2ns)']
s2 = s2.to_frame()
rmsf['CDR3a (3ns)'] = s2['CDR3a (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['CDR3a (1ns)','CDR3a (2ns)','CDR3a (3ns)']
rmsf['CDR3a'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['CDR3a SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_CDR3a = rmsf['CDR3a'].idxmax().right

print(mx_rmsf_CDR3a)

#rmsf
t, t1, t2, CDR1b, CDR1b1, CDR1b2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR1b.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR1b1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_CDR1b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR1b2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'CDR1b (1ns)': CDR1b}
df = pd.DataFrame(rmsf)
s=pd.cut(df['CDR1b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'CDR1b (2ns)': CDR1b1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['CDR1b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'CDR1b (3ns)': CDR1b2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['CDR1b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['CDR1b (2ns)'] = s1['CDR1b (2ns)']
s2 = s2.to_frame()
rmsf['CDR1b (3ns)'] = s2['CDR1b (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['CDR1b (1ns)','CDR1b (2ns)','CDR1b (3ns)']
rmsf['CDR1b'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['CDR1b SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_CDR1b = rmsf['CDR1b'].idxmax().right

print(mx_rmsf_CDR1b)


#rmsf
t, t1, t2, CDR2b, CDR2b1, CDR2b2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR2b.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR2b1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_CDR2b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR2b2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'CDR2b (1ns)': CDR2b}
df = pd.DataFrame(rmsf)
s=pd.cut(df['CDR2b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'CDR2b (2ns)': CDR2b1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['CDR2b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'CDR2b (3ns)': CDR2b2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['CDR2b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['CDR2b (2ns)'] = s1['CDR2b (2ns)']
s2 = s2.to_frame()
rmsf['CDR2b (3ns)'] = s2['CDR2b (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['CDR2b (1ns)','CDR2b (2ns)','CDR2b (3ns)']
rmsf['CDR2b'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['CDR2b SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_CDR2b = rmsf['CDR2b'].idxmax().right

print(mx_rmsf_CDR2b)

#rmsf
t, t1, t2, CDR3b, CDR3b1, CDR3b2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR3b.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR3b1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_CDR3b.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR3b2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'CDR3b (1ns)': CDR3b}
df = pd.DataFrame(rmsf)
s=pd.cut(df['CDR3b (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'CDR3b (2ns)': CDR3b1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['CDR3b (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'CDR3b (3ns)': CDR3b2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['CDR3b (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['CDR3b (2ns)'] = s1['CDR3b (2ns)']
s2 = s2.to_frame()
rmsf['CDR3b (3ns)'] = s2['CDR3b (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['CDR3b (1ns)','CDR3b (2ns)','CDR3b (3ns)']
rmsf['CDR3b'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['CDR3b SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_CDR3b = rmsf['CDR3b'].idxmax().right

print(mx_rmsf_CDR3b)

#rmsf
t, t1, t2, MHCa, MHCa1, MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            MHCa.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            MHCa1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            MHCa2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'MHCa (1ns)': MHCa}
df = pd.DataFrame(rmsf)
s=pd.cut(df['MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'MHCa (2ns)': MHCa1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'MHCa (3ns)': MHCa2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['MHCa (2ns)'] = s1['MHCa (2ns)']
s2 = s2.to_frame()
rmsf['MHCa (3ns)'] = s2['MHCa (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['MHCa (1ns)','MHCa (2ns)','MHCa (3ns)']
rmsf['MHCa'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['MHCa SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_MHCa = rmsf['MHCa'].idxmax().right

print(mx_rmsf_MHCa)

#rmsf
t, t1, t2, MHCb, MHCb1, MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            MHCb.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            MHCb1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            MHCb2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'MHCb (1ns)': MHCb}
df = pd.DataFrame(rmsf)
s=pd.cut(df['MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'MHCb (2ns)': MHCb1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'MHCb (3ns)': MHCb2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['MHCb (2ns)'] = s1['MHCb (2ns)']
s2 = s2.to_frame()
rmsf['MHCb (3ns)'] = s2['MHCb (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['MHCb (1ns)','MHCb (2ns)','MHCb (3ns)']
rmsf['MHCb'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['MHCb SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_MHCb = rmsf['MHCb'].idxmax().right

print(mx_rmsf_MHCb)

#dist
t, t1, t2, CDR3dis, CDR3dis1, CDR3dis2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/CDR3dis.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            CDR3dis.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/CDR3dis.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            CDR3dis1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/CDR3dis.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            CDR3dis2.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0.5,1.5,50).tolist()
            
dist={'CDR3dis (1ns)': CDR3dis}
df = pd.DataFrame(dist)
s=pd.cut(df['CDR3dis (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist1={'CDR3dis (2ns)': CDR3dis1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['CDR3dis (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'CDR3dis (3ns)': CDR3dis2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['CDR3dis (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['CDR3dis (2ns)'] = s1['CDR3dis (2ns)']
s2 = s2.to_frame()
dist['CDR3dis (3ns)'] = s2['CDR3dis (3ns)']

l =np.linspace(0.5,1.5,50).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['CDR3dis (1ns)','CDR3dis (2ns)','CDR3dis (3ns)']
dist['CDR3dis'] = dist[M1].astype(float).mean(axis=1)
dist['CDR3dis SEM'] = dist[M1].astype(float).sem(axis=1)

dist['ndx'] = dist.index


mx_CDR3dis = dist['CDR3dis'].idxmax().right

#instc
t, t1, t2, instc, instc1, instc2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc2.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,500,500).tolist()
            
dist={'instc (1ns)': instc}
df = pd.DataFrame(dist)
s=pd.cut(df['instc (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc (2ns)': instc1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc (3ns)': instc2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc (2ns)'] = s1['instc (2ns)']
s2 = s2.to_frame()
dist['instc (3ns)'] = s2['instc (3ns)']

l =np.linspace(0,500,500).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc (1ns)','instc (2ns)','instc (3ns)']
dist['instc'] = dist[M1].astype(float).mean(axis=1)
dist['instc SEM'] = dist[M1].astype(float).sem(axis=1)

dist['ndx'] = dist.index


mx_instc = dist['instc'].idxmax().right

#insth
t, t1, t2, insth, insth1, insth2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth2.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,50,50).tolist()
            
dist={'insth (1ns)': insth}
df = pd.DataFrame(dist)
s=pd.cut(df['insth (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth (2ns)': insth1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth (3ns)': insth2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth (2ns)'] = s1['insth (2ns)']
s2 = s2.to_frame()
dist['insth (3ns)'] = s2['insth (3ns)']

l =np.linspace(0,50,50).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth (1ns)','insth (2ns)','insth (3ns)']
dist['insth'] = dist[M1].astype(float).mean(axis=1)
dist['insth SEM'] = dist[M1].astype(float).sem(axis=1)

dist['ndx'] = dist.index


mx_insth = dist['insth'].idxmax().right

#bndx_time
t, t1, t2, bndx, bndx1, bndx2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/pullxcf.xvg") as f:
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
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/pullxcf.xvg") as f:
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
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            bndx2.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(5,15,100).tolist()
            
dist={'bndx (1ns)': bndx}
df = pd.DataFrame(dist)
s=pd.cut(df['bndx (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist1={'bndx (2ns)': bndx1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['bndx (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'bndx (3ns)': bndx2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['bndx (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['bndx (2ns)'] = s1['bndx (2ns)']
s2 = s2.to_frame()
dist['bndx (3ns)'] = s2['bndx (3ns)']

l =np.linspace(5,15,100).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['bndx (1ns)','bndx (2ns)','bndx (3ns)']
dist['bndx'] = dist[M1].astype(float).mean(axis=1)
dist['bndx SEM'] = dist[M1].astype(float).sem(axis=1)

dist['ndx'] = dist.index

mx_bndx = dist['bndx'].idxmax().right

time1 = len(t)
time2 = len(t1)
time3 = len(t2)

time = (time1 + time2 + time3)/30


#RXN Distance

df = pd.read_excel (r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/iex.xlsx')
df2 = pd.read_excel (r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/iex.xlsx')
df3 = pd.read_excel (r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/iex.xlsx')

df = df.append(df2, ignore_index=True)
df = df.append(df3, ignore_index=True)

bins=np.linspace(5.884,14.612,200).tolist()
MHCc = df.groupby(pd.cut(df['COM Distance'], bins=bins))['MHC-Coulombic'].agg(['mean','sem','size'])
MHClj = df.groupby(pd.cut(df['COM Distance'], bins=bins))['MHC-LennardJones'].agg(['mean','sem','size'])
pepc = df.groupby(pd.cut(df['COM Distance'], bins=bins))['Peptide_Coulombic'].agg(['mean','sem','size'])
peplj = df.groupby(pd.cut(df['COM Distance'], bins=bins))['Peptide_LennardJones'].agg(['mean','sem','size'])


zero_MHCc = MHCc['mean'].idxmax().right
zero_MHClj = MHClj['mean'].idxmax().right
zero_pepc = pepc['mean'].idxmax().right
zero_peplj = peplj['mean'].idxmax().right

zero = max(zero_MHCc,zero_MHClj,zero_pepc,zero_peplj)

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2= [],[],[],[],[],[]
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts.log") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            d.append(str(cols[0]))
            a.append(str(cols[1]))

# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts.log") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            d1.append(str(cols[0]))
            a1.append(str(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts.log") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            d2.append(str(cols[0]))
            a2.append(str(cols[1]))



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts = (con1 + con2 + con3)/3


# Total H-Bonds

hmap1 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/hmap.csv', index_col = 0, header=None)
hmap2 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/hmap.csv', index_col = 0, header=None)
hmap3 = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/hmap.csv', index_col = 0, header=None)

hb1, t1 = hmap1.shape
hb2, t2 = hmap2.shape
hb3, t3 =  hmap3.shape

hb = (hb1 + hb2 + hb3)/3

#insth
t, t1, t2, insth_CDR1a_MHCa, insth_CDR1a_MHCa1, insth_CDR1a_MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR1a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR1a_MHCa.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR1a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR1a_MHCa1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR1a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR1a_MHCa1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR1a_MHCa (1ns)': insth_CDR1a_MHCa}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR1a_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR1a_MHCa (2ns)': insth_CDR1a_MHCa1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR1a_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR1a_MHCa (3ns)': insth_CDR1a_MHCa2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR1a_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR1a_MHCa (2ns)'] = s1['insth_CDR1a_MHCa (2ns)']
s2 = s2.to_frame()
dist['insth_CDR1a_MHCa (3ns)'] = s2['insth_CDR1a_MHCa (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR1a_MHCa (1ns)','insth_CDR1a_MHCa (2ns)','insth_CDR1a_MHCa (3ns)']
dist['insth_CDR1a_MHCa'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR1a_MHCa SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR1a_MHCa = dist['insth_CDR1a_MHCa'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR1a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR1a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')                
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR1a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')


con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR1a_MHCa = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR1a_pep, insth_CDR1a_pep1, insth_CDR1a_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR1a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR1a_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR1a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR1a_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR1a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR1a_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR1a_pep (1ns)': insth_CDR1a_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR1a_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR1a_pep (2ns)': insth_CDR1a_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR1a_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR1a_pep (3ns)': insth_CDR1a_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR1a_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR1a_pep (2ns)'] = s1['insth_CDR1a_pep (2ns)']
s2 = s2.to_frame()
dist['insth_CDR1a_pep (3ns)'] = s2['insth_CDR1a_pep (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR1a_pep (1ns)','insth_CDR1a_pep (2ns)','insth_CDR1a_pep (3ns)']
dist['insth_CDR1a_pep'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR1a_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR1a_pep = dist['insth_CDR1a_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR1a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR1a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR1a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR1a_pep = (con1 + con2 + con3)/3

#insth
t, t1, t2, insth_CDR2a_MHCa, insth_CDR2a_MHCa1, insth_CDR2a_MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR2a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR2a_MHCa.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR2a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR2a_MHCa1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR2a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR2a_MHCa1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR2a_MHCa (1ns)': insth_CDR2a_MHCa}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR2a_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR2a_MHCa (2ns)': insth_CDR2a_MHCa1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR2a_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR2a_MHCa (3ns)': insth_CDR2a_MHCa2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR2a_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR2a_MHCa (2ns)'] = s1['insth_CDR2a_MHCa (2ns)']
s2 = s2.to_frame()
dist['insth_CDR2a_MHCa (3ns)'] = s2['insth_CDR2a_MHCa (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR2a_MHCa (1ns)','insth_CDR2a_MHCa (2ns)','insth_CDR2a_MHCa (3ns)']
dist['insth_CDR2a_MHCa'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR2a_MHCa SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR2a_MHCa = dist['insth_CDR2a_MHCa'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR2a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR2a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR2a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR2a_MHCa = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR2a_pep, insth_CDR2a_pep1, insth_CDR2a_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR2a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR2a_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR2a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR2a_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR2a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR2a_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR2a_pep (1ns)': insth_CDR2a_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR2a_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR2a_pep (2ns)': insth_CDR2a_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR2a_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR2a_pep (3ns)': insth_CDR2a_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR2a_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR2a_pep (2ns)'] = s1['insth_CDR2a_pep (2ns)']
s2 = s2.to_frame()
dist['insth_CDR2a_pep (3ns)'] = s2['insth_CDR2a_pep (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR2a_pep (1ns)','insth_CDR2a_pep (2ns)','insth_CDR2a_pep (3ns)']
dist['insth_CDR2a_pep'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR2a_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR2a_pep = dist['insth_CDR2a_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR2a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')
    
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR2a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
                
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR2a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')                



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR2a_pep = (con1 + con2 + con3)/3

#insth
t, t1, t2, insth_CDR3a_MHCa, insth_CDR3a_MHCa1, insth_CDR3a_MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR3a_MHCa.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t.append(float(cols[0]))
                insth_CDR3a_MHCa.append(float(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR3a_MHCa.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t1.append(float(cols[0]))
                insth_CDR3a_MHCa1.append(float(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR3a_MHCa.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t2.append(float(cols[0]))
                insth_CDR3a_MHCa1.append(float(cols[1]))
except Exception:
    print('not found')

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR3a_MHCa (1ns)': insth_CDR3a_MHCa}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR3a_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR3a_MHCa (2ns)': insth_CDR3a_MHCa1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR3a_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR3a_MHCa (3ns)': insth_CDR3a_MHCa2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR3a_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR3a_MHCa (2ns)'] = s1['insth_CDR3a_MHCa (2ns)']
s2 = s2.to_frame()
dist['insth_CDR3a_MHCa (3ns)'] = s2['insth_CDR3a_MHCa (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR3a_MHCa (1ns)','insth_CDR3a_MHCa (2ns)','insth_CDR3a_MHCa (3ns)']
dist['insth_CDR3a_MHCa'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR3a_MHCa SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index



mx_insth_CDR3a_MHCa = dist['insth_CDR3a_MHCa'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR3a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR3a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR3a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR3a_MHCa = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR3a_pep, insth_CDR3a_pep1, insth_CDR3a_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR3a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR3a_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR3a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR3a_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR3a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR3a_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR3a_pep (1ns)': insth_CDR3a_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR3a_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR3a_pep (2ns)': insth_CDR3a_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR3a_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR3a_pep (3ns)': insth_CDR3a_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR3a_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR3a_pep (2ns)'] = s1['insth_CDR3a_pep (2ns)']
s2 = s2.to_frame()
dist['insth_CDR3a_pep (3ns)'] = s2['insth_CDR3a_pep (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR3a_pep (1ns)','insth_CDR3a_pep (2ns)','insth_CDR3a_pep (3ns)']
dist['insth_CDR3a_pep'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR3a_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index

mx_insth_CDR3a_pep = dist['insth_CDR3a_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR3a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR3a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR3a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR3a_pep = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR1b_MHCb, insth_CDR1b_MHCb1, insth_CDR1b_MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR1b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR1b_MHCb.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR1b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR1b_MHCb1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR1b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR1b_MHCb1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR1b_MHCb (1ns)': insth_CDR1b_MHCb}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR1b_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR1b_MHCb (2ns)': insth_CDR1b_MHCb1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR1b_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR1b_MHCb (3ns)': insth_CDR1b_MHCb2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR1b_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR1b_MHCb (2ns)'] = s1['insth_CDR1b_MHCb (2ns)']
s2 = s2.to_frame()
dist['insth_CDR1b_MHCb (3ns)'] = s2['insth_CDR1b_MHCb (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR1b_MHCb (1ns)','insth_CDR1b_MHCb (2ns)','insth_CDR1b_MHCb (3ns)']
dist['insth_CDR1b_MHCb'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR1b_MHCb SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR1b_MHCb = dist['insth_CDR1b_MHCb'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR1b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR1b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')                
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR1b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')


con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR1b_MHCb = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR1b_pep, insth_CDR1b_pep1, insth_CDR1b_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR1b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR1b_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR1b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR1b_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR1b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR1b_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR1b_pep (1ns)': insth_CDR1b_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR1b_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR1b_pep (2ns)': insth_CDR1b_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR1b_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR1b_pep (3ns)': insth_CDR1b_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR1b_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR1b_pep (2ns)'] = s1['insth_CDR1b_pep (2ns)']
s2 = s2.to_frame()
dist['insth_CDR1b_pep (3ns)'] = s2['insth_CDR1b_pep (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR1b_pep (1ns)','insth_CDR1b_pep (2ns)','insth_CDR1b_pep (3ns)']
dist['insth_CDR1b_pep'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR1b_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR1b_pep = dist['insth_CDR1b_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR1b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR1b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR1b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR1b_pep = (con1 + con2 + con3)/3

#insth
t, t1, t2, insth_CDR2b_MHCb, insth_CDR2b_MHCb1, insth_CDR2b_MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR2b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR2b_MHCb.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR2b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR2b_MHCb1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR2b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR2b_MHCb1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR2b_MHCb (1ns)': insth_CDR2b_MHCb}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR2b_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR2b_MHCb (2ns)': insth_CDR2b_MHCb1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR2b_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR2b_MHCb (3ns)': insth_CDR2b_MHCb2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR2b_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR2b_MHCb (2ns)'] = s1['insth_CDR2b_MHCb (2ns)']
s2 = s2.to_frame()
dist['insth_CDR2b_MHCb (3ns)'] = s2['insth_CDR2b_MHCb (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR2b_MHCb (1ns)','insth_CDR2b_MHCb (2ns)','insth_CDR2b_MHCb (3ns)']
dist['insth_CDR2b_MHCb'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR2b_MHCb SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR2b_MHCb = dist['insth_CDR2b_MHCb'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR2b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR2b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR2b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR2b_MHCb = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR2b_pep, insth_CDR2b_pep1, insth_CDR2b_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR2b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR2b_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR2b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR2b_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR2b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR2b_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR2b_pep (1ns)': insth_CDR2b_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR2b_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR2b_pep (2ns)': insth_CDR2b_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR2b_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR2b_pep (3ns)': insth_CDR2b_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR2b_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR2b_pep (2ns)'] = s1['insth_CDR2b_pep (2ns)']
s2 = s2.to_frame()
dist['insth_CDR2b_pep (3ns)'] = s2['insth_CDR2b_pep (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR2b_pep (1ns)','insth_CDR2b_pep (2ns)','insth_CDR2b_pep (3ns)']
dist['insth_CDR2b_pep'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR2b_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_insth_CDR2b_pep = dist['insth_CDR2b_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR2b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')
    
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR2b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
                
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR2b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')                



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR2b_pep = (con1 + con2 + con3)/3

#insth
t, t1, t2, insth_CDR3b_MHCb, insth_CDR3b_MHCb1, insth_CDR3b_MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR3b_MHCb.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t.append(float(cols[0]))
                insth_CDR3b_MHCb.append(float(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR3b_MHCb.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t1.append(float(cols[0]))
                insth_CDR3b_MHCb1.append(float(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR3b_MHCb.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t2.append(float(cols[0]))
                insth_CDR3b_MHCb1.append(float(cols[1]))
except Exception:
    print('not found')

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR3b_MHCb (1ns)': insth_CDR3b_MHCb}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR3b_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR3b_MHCb (2ns)': insth_CDR3b_MHCb1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR3b_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR3b_MHCb (3ns)': insth_CDR3b_MHCb2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR3b_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR3b_MHCb (2ns)'] = s1['insth_CDR3b_MHCb (2ns)']
s2 = s2.to_frame()
dist['insth_CDR3b_MHCb (3ns)'] = s2['insth_CDR3b_MHCb (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR3b_MHCb (1ns)','insth_CDR3b_MHCb (2ns)','insth_CDR3b_MHCb (3ns)']
dist['insth_CDR3b_MHCb'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR3b_MHCb SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index



mx_insth_CDR3b_MHCb = dist['insth_CDR3b_MHCb'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR3b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR3b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR3b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR3b_MHCb = (con1 + con2 + con3)/3


#insth
t, t1, t2, insth_CDR3b_pep, insth_CDR3b_pep1, insth_CDR3b_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/num_interface_CDR3b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            insth_CDR3b_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/num_interface_CDR3b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            insth_CDR3b_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/num_interface_CDR3b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            insth_CDR3b_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,19,39).tolist()
            
dist={'insth_CDR3b_pep (1ns)': insth_CDR3b_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['insth_CDR3b_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'insth_CDR3b_pep (2ns)': insth_CDR3b_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['insth_CDR3b_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'insth_CDR3b_pep (3ns)': insth_CDR3b_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['insth_CDR3b_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['insth_CDR3b_pep (2ns)'] = s1['insth_CDR3b_pep (2ns)']
s2 = s2.to_frame()
dist['insth_CDR3b_pep (3ns)'] = s2['insth_CDR3b_pep (3ns)']

l =np.linspace(0,19,39).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['insth_CDR3b_pep (1ns)','insth_CDR3b_pep (2ns)','insth_CDR3b_pep (3ns)']
dist['insth_CDR3b_pep'] = dist[M1].astype(float).mean(axis=1)
dist['insth_CDR3b_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index

mx_insth_CDR3b_pep = dist['insth_CDR3b_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/h_CDR3b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
                h.append(str(cols[2]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/h_CDR3b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
                h1.append(str(cols[2]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/h_CDR3b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
                h2.append(str(cols[2]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

h_CDR3b_pep = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR1a_MHCa, instc_CDR1a_MHCa1, instc_CDR1a_MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR1a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR1a_MHCa.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR1a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR1a_MHCa1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR1a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR1a_MHCa1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR1a_MHCa (1ns)': instc_CDR1a_MHCa}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR1a_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR1a_MHCa (2ns)': instc_CDR1a_MHCa1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR1a_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR1a_MHCa (3ns)': instc_CDR1a_MHCa2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR1a_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR1a_MHCa (2ns)'] = s1['instc_CDR1a_MHCa (2ns)']
s2 = s2.to_frame()
dist['instc_CDR1a_MHCa (3ns)'] = s2['instc_CDR1a_MHCa (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR1a_MHCa (1ns)','instc_CDR1a_MHCa (2ns)','instc_CDR1a_MHCa (3ns)']
dist['instc_CDR1a_MHCa'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR1a_MHCa SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR1a_MHCa = dist['instc_CDR1a_MHCa'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR1a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR1a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')                
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR1a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')


con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR1a_MHCa = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR1a_pep, instc_CDR1a_pep1, instc_CDR1a_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR1a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR1a_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR1a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR1a_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR1a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR1a_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR1a_pep (1ns)': instc_CDR1a_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR1a_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR1a_pep (2ns)': instc_CDR1a_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR1a_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR1a_pep (3ns)': instc_CDR1a_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR1a_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR1a_pep (2ns)'] = s1['instc_CDR1a_pep (2ns)']
s2 = s2.to_frame()
dist['instc_CDR1a_pep (3ns)'] = s2['instc_CDR1a_pep (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR1a_pep (1ns)','instc_CDR1a_pep (2ns)','instc_CDR1a_pep (3ns)']
dist['instc_CDR1a_pep'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR1a_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR1a_pep = dist['instc_CDR1a_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR1a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR1a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR1a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR1a_pep = (con1 + con2 + con3)/3

#instc
t, t1, t2, instc_CDR2a_MHCa, instc_CDR2a_MHCa1, instc_CDR2a_MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR2a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR2a_MHCa.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR2a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR2a_MHCa1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR2a_MHCa.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR2a_MHCa1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR2a_MHCa (1ns)': instc_CDR2a_MHCa}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR2a_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR2a_MHCa (2ns)': instc_CDR2a_MHCa1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR2a_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR2a_MHCa (3ns)': instc_CDR2a_MHCa2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR2a_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR2a_MHCa (2ns)'] = s1['instc_CDR2a_MHCa (2ns)']
s2 = s2.to_frame()
dist['instc_CDR2a_MHCa (3ns)'] = s2['instc_CDR2a_MHCa (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR2a_MHCa (1ns)','instc_CDR2a_MHCa (2ns)','instc_CDR2a_MHCa (3ns)']
dist['instc_CDR2a_MHCa'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR2a_MHCa SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR2a_MHCa = dist['instc_CDR2a_MHCa'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR2a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR2a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR2a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR2a_MHCa = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR2a_pep, instc_CDR2a_pep1, instc_CDR2a_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR2a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR2a_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR2a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR2a_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR2a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR2a_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR2a_pep (1ns)': instc_CDR2a_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR2a_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR2a_pep (2ns)': instc_CDR2a_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR2a_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR2a_pep (3ns)': instc_CDR2a_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR2a_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR2a_pep (2ns)'] = s1['instc_CDR2a_pep (2ns)']
s2 = s2.to_frame()
dist['instc_CDR2a_pep (3ns)'] = s2['instc_CDR2a_pep (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR2a_pep (1ns)','instc_CDR2a_pep (2ns)','instc_CDR2a_pep (3ns)']
dist['instc_CDR2a_pep'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR2a_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR2a_pep = dist['instc_CDR2a_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR2a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')
    
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR2a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
                
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR2a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')                



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR2a_pep = (con1 + con2 + con3)/3

#instc
t, t1, t2, instc_CDR3a_MHCa, instc_CDR3a_MHCa1, instc_CDR3a_MHCa2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR3a_MHCa.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t.append(float(cols[0]))
                instc_CDR3a_MHCa.append(float(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR3a_MHCa.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t1.append(float(cols[0]))
                instc_CDR3a_MHCa1.append(float(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR3a_MHCa.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t2.append(float(cols[0]))
                instc_CDR3a_MHCa1.append(float(cols[1]))
except Exception:
    print('not found')

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR3a_MHCa (1ns)': instc_CDR3a_MHCa}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR3a_MHCa (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR3a_MHCa (2ns)': instc_CDR3a_MHCa1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR3a_MHCa (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR3a_MHCa (3ns)': instc_CDR3a_MHCa2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR3a_MHCa (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR3a_MHCa (2ns)'] = s1['instc_CDR3a_MHCa (2ns)']
s2 = s2.to_frame()
dist['instc_CDR3a_MHCa (3ns)'] = s2['instc_CDR3a_MHCa (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR3a_MHCa (1ns)','instc_CDR3a_MHCa (2ns)','instc_CDR3a_MHCa (3ns)']
dist['instc_CDR3a_MHCa'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR3a_MHCa SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index



mx_instc_CDR3a_MHCa = dist['instc_CDR3a_MHCa'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR3a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR3a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR3a_MHCa.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR3a_MHCa = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR3a_pep, instc_CDR3a_pep1, instc_CDR3a_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR3a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR3a_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR3a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR3a_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR3a_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR3a_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR3a_pep (1ns)': instc_CDR3a_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR3a_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR3a_pep (2ns)': instc_CDR3a_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR3a_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR3a_pep (3ns)': instc_CDR3a_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR3a_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR3a_pep (2ns)'] = s1['instc_CDR3a_pep (2ns)']
s2 = s2.to_frame()
dist['instc_CDR3a_pep (3ns)'] = s2['instc_CDR3a_pep (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR3a_pep (1ns)','instc_CDR3a_pep (2ns)','instc_CDR3a_pep (3ns)']
dist['instc_CDR3a_pep'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR3a_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index

mx_instc_CDR3a_pep = dist['instc_CDR3a_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR3a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR3a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR3a_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR3a_pep = (con1 + con2 + con3)/3

#instc
t, t1, t2, instc_CDR1b_MHCb, instc_CDR1b_MHCb1, instc_CDR1b_MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR1b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR1b_MHCb.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR1b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR1b_MHCb1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR1b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR1b_MHCb1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR1b_MHCb (1ns)': instc_CDR1b_MHCb}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR1b_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR1b_MHCb (2ns)': instc_CDR1b_MHCb1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR1b_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR1b_MHCb (3ns)': instc_CDR1b_MHCb2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR1b_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR1b_MHCb (2ns)'] = s1['instc_CDR1b_MHCb (2ns)']
s2 = s2.to_frame()
dist['instc_CDR1b_MHCb (3ns)'] = s2['instc_CDR1b_MHCb (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR1b_MHCb (1ns)','instc_CDR1b_MHCb (2ns)','instc_CDR1b_MHCb (3ns)']
dist['instc_CDR1b_MHCb'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR1b_MHCb SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR1b_MHCb = dist['instc_CDR1b_MHCb'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR1b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR1b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')                
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR1b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')


con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR1b_MHCb = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR1b_pep, instc_CDR1b_pep1, instc_CDR1b_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR1b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR1b_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR1b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR1b_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR1b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR1b_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR1b_pep (1ns)': instc_CDR1b_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR1b_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR1b_pep (2ns)': instc_CDR1b_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR1b_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR1b_pep (3ns)': instc_CDR1b_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR1b_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR1b_pep (2ns)'] = s1['instc_CDR1b_pep (2ns)']
s2 = s2.to_frame()
dist['instc_CDR1b_pep (3ns)'] = s2['instc_CDR1b_pep (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR1b_pep (1ns)','instc_CDR1b_pep (2ns)','instc_CDR1b_pep (3ns)']
dist['instc_CDR1b_pep'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR1b_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR1b_pep = dist['instc_CDR1b_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR1b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR1b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR1b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR1b_pep = (con1 + con2 + con3)/3

#instc
t, t1, t2, instc_CDR2b_MHCb, instc_CDR2b_MHCb1, instc_CDR2b_MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR2b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR2b_MHCb.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR2b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR2b_MHCb1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR2b_MHCb.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR2b_MHCb1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR2b_MHCb (1ns)': instc_CDR2b_MHCb}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR2b_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR2b_MHCb (2ns)': instc_CDR2b_MHCb1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR2b_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR2b_MHCb (3ns)': instc_CDR2b_MHCb2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR2b_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR2b_MHCb (2ns)'] = s1['instc_CDR2b_MHCb (2ns)']
s2 = s2.to_frame()
dist['instc_CDR2b_MHCb (3ns)'] = s2['instc_CDR2b_MHCb (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR2b_MHCb (1ns)','instc_CDR2b_MHCb (2ns)','instc_CDR2b_MHCb (3ns)']
dist['instc_CDR2b_MHCb'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR2b_MHCb SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR2b_MHCb = dist['instc_CDR2b_MHCb'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR2b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR2b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR2b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR2b_MHCb = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR2b_pep, instc_CDR2b_pep1, instc_CDR2b_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR2b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR2b_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR2b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR2b_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR2b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR2b_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR2b_pep (1ns)': instc_CDR2b_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR2b_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR2b_pep (2ns)': instc_CDR2b_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR2b_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR2b_pep (3ns)': instc_CDR2b_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR2b_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR2b_pep (2ns)'] = s1['instc_CDR2b_pep (2ns)']
s2 = s2.to_frame()
dist['instc_CDR2b_pep (3ns)'] = s2['instc_CDR2b_pep (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR2b_pep (1ns)','instc_CDR2b_pep (2ns)','instc_CDR2b_pep (3ns)']
dist['instc_CDR2b_pep'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR2b_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index


mx_instc_CDR2b_pep = dist['instc_CDR2b_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR2b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')
    
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR2b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
                
# 3ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR2b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')                



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR2b_pep = (con1 + con2 + con3)/3

#instc
t, t1, t2, instc_CDR3b_MHCb, instc_CDR3b_MHCb1, instc_CDR3b_MHCb2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns
try: 
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR3b_MHCb.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t.append(float(cols[0]))
                instc_CDR3b_MHCb.append(float(cols[1]))
except Exception:
    print('not found')
# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR3b_MHCb.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t1.append(float(cols[0]))
                instc_CDR3b_MHCb1.append(float(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR3b_MHCb.xvg") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 3:
                t2.append(float(cols[0]))
                instc_CDR3b_MHCb1.append(float(cols[1]))
except Exception:
    print('not found')

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR3b_MHCb (1ns)': instc_CDR3b_MHCb}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR3b_MHCb (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR3b_MHCb (2ns)': instc_CDR3b_MHCb1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR3b_MHCb (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR3b_MHCb (3ns)': instc_CDR3b_MHCb2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR3b_MHCb (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR3b_MHCb (2ns)'] = s1['instc_CDR3b_MHCb (2ns)']
s2 = s2.to_frame()
dist['instc_CDR3b_MHCb (3ns)'] = s2['instc_CDR3b_MHCb (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR3b_MHCb (1ns)','instc_CDR3b_MHCb (2ns)','instc_CDR3b_MHCb (3ns)']
dist['instc_CDR3b_MHCb'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR3b_MHCb SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index



mx_instc_CDR3b_MHCb = dist['instc_CDR3b_MHCb'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR3b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR3b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR3b_MHCb.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR3b_MHCb = (con1 + con2 + con3)/3


#instc
t, t1, t2, instc_CDR3b_pep, instc_CDR3b_pep1, instc_CDR3b_pep2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/numc_interface_CDR3b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t.append(float(cols[0]))
            instc_CDR3b_pep.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/numc_interface_CDR3b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t1.append(float(cols[0]))
            instc_CDR3b_pep1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/numc_interface_CDR3b_pep.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            t2.append(float(cols[0]))
            instc_CDR3b_pep1.append(float(cols[1]))

#Bin dist  
cut_bins=np.linspace(0,49,99).tolist()
            
dist={'instc_CDR3b_pep (1ns)': instc_CDR3b_pep}
df = pd.DataFrame(dist)
s=pd.cut(df['instc_CDR3b_pep (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)
print(s)
dist1={'instc_CDR3b_pep (2ns)': instc_CDR3b_pep1}
df1 = pd.DataFrame(dist1)
s1=pd.cut(df1['instc_CDR3b_pep (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist2={'instc_CDR3b_pep (3ns)': instc_CDR3b_pep2}
df2 = pd.DataFrame(dist2)
s2=pd.cut(df2['instc_CDR3b_pep (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

dist = s.to_frame()
s1 = s1.to_frame()
dist['instc_CDR3b_pep (2ns)'] = s1['instc_CDR3b_pep (2ns)']
s2 = s2.to_frame()
dist['instc_CDR3b_pep (3ns)'] = s2['instc_CDR3b_pep (3ns)']

l =np.linspace(0,49,99).tolist()
x = pd.Series(l[:-1])
dist['dist']= x.values

M1 = ['instc_CDR3b_pep (1ns)','instc_CDR3b_pep (2ns)','instc_CDR3b_pep (3ns)']
dist['instc_CDR3b_pep'] = dist[M1].astype(float).mean(axis=1)
dist['instc_CDR3b_pep SEM'] = dist[M1].astype(float).sem(axis=1)
dist = dist.fillna(0)
dist['ndx'] = dist.index

mx_instc_CDR3b_pep = dist['instc_CDR3b_pep'].idxmax().right

#total contacts
##Data Acquisition
# 1ns 

d, d1, d2, a, a1, a2, h, h1, h2 = [],[],[],[],[],[],[],[],[]
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/contacts_CDR3b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d.append(str(cols[0]))
                a.append(str(cols[1]))
except Exception:
    print('not found')

# 2ns
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/contacts_CDR3b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d1.append(str(cols[0]))
                a1.append(str(cols[1]))
except Exception:
    print('not found')
# 3ns 
try:
    with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/contacts_CDR3b_pep.log") as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                continue
            cols = line.split()
    
            if len(cols) == 2:
                d2.append(str(cols[0]))
                a2.append(str(cols[1]))
except Exception:
    print('not found')



con1 = len(d)
con2 = len(d1)
con3 = len(d2)

contacts_CDR3b_pep = (con1 + con2 + con3)/3

#gyr
t, t1, t2, gyr_TCR, gyr_TCR1, gyr_TCR2, gyrx_TCR, gyrx_TCR1, gyrx_TCR2, gyry_TCR, gyry_TCR1, gyry_TCR2, gyrz_TCR, gyrz_TCR1, gyrz_TCR2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_TCR.append(float(cols[1]))
            gyrx_TCR.append(float(cols[2]))
            gyry_TCR.append(float(cols[3]))
            gyrz_TCR.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_TCR1.append(float(cols[1]))
            gyrx_TCR1.append(float(cols[2]))
            gyry_TCR1.append(float(cols[3]))
            gyrz_TCR1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_TCR2.append(float(cols[1]))
            gyrx_TCR2.append(float(cols[2]))
            gyry_TCR2.append(float(cols[3]))
            gyrz_TCR2.append(float(cols[4]))

#Bin gyr  

cut_bins=np.linspace(0,3,300).tolist()
            
gyr={'gyr_TCR (1ns)': gyr_TCR}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_TCR (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_TCR (2ns)': gyr_TCR1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_TCR (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_TCR (3ns)': gyr_TCR2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_TCR (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_TCR (2ns)'] = s1['gyr_TCR (2ns)']
s2 = s2.to_frame()
gyr['gyr_TCR (3ns)'] = s2['gyr_TCR (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_TCR (1ns)','gyr_TCR (2ns)','gyr_TCR (3ns)']
gyr['gyr_TCR'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_TCR SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_TCR = gyr['gyr_TCR'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,3,300).tolist()
            
gyrx={'gyrx_TCR (1ns)': gyrx_TCR}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_TCR (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_TCR (2ns)': gyrx_TCR1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_TCR (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_TCR (3ns)': gyrx_TCR2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_TCR (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_TCR (2ns)'] = s1['gyrx_TCR (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_TCR (3ns)'] = s2['gyrx_TCR (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_TCR (1ns)','gyrx_TCR (2ns)','gyrx_TCR (3ns)']
gyrx['gyrx_TCR'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_TCR SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_TCR = gyrx['gyrx_TCR'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,3,300).tolist()
            
gyry={'gyry_TCR (1ns)': gyry_TCR}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_TCR (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_TCR (2ns)': gyry_TCR1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_TCR (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_TCR (3ns)': gyry_TCR2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_TCR (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_TCR (2ns)'] = s1['gyry_TCR (2ns)']
s2 = s2.to_frame()
gyry['gyry_TCR (3ns)'] = s2['gyry_TCR (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_TCR (1ns)','gyry_TCR (2ns)','gyry_TCR (3ns)']
gyry['gyry_TCR'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_TCR SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_TCR = gyry['gyry_TCR'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,3,300).tolist()
            
gyrz={'gyrz_TCR (1ns)': gyrz_TCR}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_TCR (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_TCR (2ns)': gyrz_TCR1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_TCR (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_TCR (3ns)': gyrz_TCR2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_TCR (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_TCR (2ns)'] = s1['gyrz_TCR (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_TCR (3ns)'] = s2['gyrz_TCR (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_TCR (1ns)','gyrz_TCR (2ns)','gyrz_TCR (3ns)']
gyrz['gyrz_TCR'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_TCR SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_TCR = gyrz['gyrz_TCR'].idxmax().right

#gyr
t, t1, t2, gyr_pMHC, gyr_pMHC1, gyr_pMHC2, gyrx_pMHC, gyrx_pMHC1, gyrx_pMHC2, gyry_pMHC, gyry_pMHC1, gyry_pMHC2, gyrz_pMHC, gyrz_pMHC1, gyrz_pMHC2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/gyr_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t.append(float(cols[0]))
            gyr_pMHC.append(float(cols[1]))
            gyrx_pMHC.append(float(cols[2]))
            gyry_pMHC.append(float(cols[3]))
            gyrz_pMHC.append(float(cols[4]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/gyr_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t1.append(float(cols[0]))
            gyr_pMHC1.append(float(cols[1]))
            gyrx_pMHC1.append(float(cols[2]))
            gyry_pMHC1.append(float(cols[3]))
            gyrz_pMHC1.append(float(cols[4]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/gyr_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 5:
            t2.append(float(cols[0]))
            gyr_pMHC2.append(float(cols[1]))
            gyrx_pMHC2.append(float(cols[2]))
            gyry_pMHC2.append(float(cols[3]))
            gyrz_pMHC2.append(float(cols[4]))

#Bin gyr  

cut_bins=np.linspace(0,3,300).tolist()
            
gyr={'gyr_pMHC (1ns)': gyr_pMHC}
df = pd.DataFrame(gyr)
s=pd.cut(df['gyr_pMHC (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr1={'gyr_pMHC (2ns)': gyr_pMHC1}
df1 = pd.DataFrame(gyr1)
s1=pd.cut(df1['gyr_pMHC (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr2={'gyr_pMHC (3ns)': gyr_pMHC2}
df2 = pd.DataFrame(gyr2)
s2=pd.cut(df2['gyr_pMHC (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyr = s.to_frame()
s1 = s1.to_frame()
gyr['gyr_pMHC (2ns)'] = s1['gyr_pMHC (2ns)']
s2 = s2.to_frame()
gyr['gyr_pMHC (3ns)'] = s2['gyr_pMHC (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyr['gyr']= x.values

M1 = ['gyr_pMHC (1ns)','gyr_pMHC (2ns)','gyr_pMHC (3ns)']
gyr['gyr_pMHC'] = gyr[M1].astype(float).mean(axis=1)
gyr['gyr_pMHC SEM'] = gyr[M1].astype(float).sem(axis=1)

gyr['ndx'] = gyr.index
#print(gyr)
mx_gyr_pMHC = gyr['gyr_pMHC'].idxmax().right


#Bin gyrx  
cut_bins=np.linspace(0,3,300).tolist()
            
gyrx={'gyrx_pMHC (1ns)': gyrx_pMHC}
df = pd.DataFrame(gyrx)
s=pd.cut(df['gyrx_pMHC (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx1={'gyrx_pMHC (2ns)': gyrx_pMHC1}
df1 = pd.DataFrame(gyrx1)
s1=pd.cut(df1['gyrx_pMHC (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx2={'gyrx_pMHC (3ns)': gyrx_pMHC2}
df2 = pd.DataFrame(gyrx2)
s2=pd.cut(df2['gyrx_pMHC (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrx = s.to_frame()
s1 = s1.to_frame()
gyrx['gyrx_pMHC (2ns)'] = s1['gyrx_pMHC (2ns)']
s2 = s2.to_frame()
gyrx['gyrx_pMHC (3ns)'] = s2['gyrx_pMHC (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyrx['gyrx']= x.values

M1 = ['gyrx_pMHC (1ns)','gyrx_pMHC (2ns)','gyrx_pMHC (3ns)']
gyrx['gyrx_pMHC'] = gyrx[M1].astype(float).mean(axis=1)
gyrx['gyrx_pMHC SEM'] = gyrx[M1].astype(float).sem(axis=1)

gyrx['ndx'] = gyrx.index

mx_gyrx_pMHC = gyrx['gyrx_pMHC'].idxmax().right

#Bin gyry  
cut_bins=np.linspace(0,3,300).tolist()
            
gyry={'gyry_pMHC (1ns)': gyry_pMHC}
df = pd.DataFrame(gyry)
s=pd.cut(df['gyry_pMHC (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry1={'gyry_pMHC (2ns)': gyry_pMHC1}
df1 = pd.DataFrame(gyry1)
s1=pd.cut(df1['gyry_pMHC (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry2={'gyry_pMHC (3ns)': gyry_pMHC2}
df2 = pd.DataFrame(gyry2)
s2=pd.cut(df2['gyry_pMHC (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyry = s.to_frame()
s1 = s1.to_frame()
gyry['gyry_pMHC (2ns)'] = s1['gyry_pMHC (2ns)']
s2 = s2.to_frame()
gyry['gyry_pMHC (3ns)'] = s2['gyry_pMHC (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyry['gyry']= x.values

M1 = ['gyry_pMHC (1ns)','gyry_pMHC (2ns)','gyry_pMHC (3ns)']
gyry['gyry_pMHC'] = gyry[M1].astype(float).mean(axis=1)
gyry['gyry_pMHC SEM'] = gyry[M1].astype(float).sem(axis=1)

gyry['ndx'] = gyry.index

mx_gyry_pMHC = gyry['gyry_pMHC'].idxmax().right

#Bin gyrz  
cut_bins=np.linspace(0,3,300).tolist()
            
gyrz={'gyrz_pMHC (1ns)': gyrz_pMHC}
df = pd.DataFrame(gyrz)
s=pd.cut(df['gyrz_pMHC (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz1={'gyrz_pMHC (2ns)': gyrz_pMHC1}
df1 = pd.DataFrame(gyrz1)
s1=pd.cut(df1['gyrz_pMHC (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz2={'gyrz_pMHC (3ns)': gyrz_pMHC2}
df2 = pd.DataFrame(gyrz2)
s2=pd.cut(df2['gyrz_pMHC (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

gyrz = s.to_frame()
s1 = s1.to_frame()
gyrz['gyrz_pMHC (2ns)'] = s1['gyrz_pMHC (2ns)']
s2 = s2.to_frame()
gyrz['gyrz_pMHC (3ns)'] = s2['gyrz_pMHC (3ns)']

l =np.linspace(0,3,300).tolist()
x = pd.Series(l[:-1])
gyrz['gyrz']= x.values

M1 = ['gyrz_pMHC (1ns)','gyrz_pMHC (2ns)','gyrz_pMHC (3ns)']
gyrz['gyrz_pMHC'] = gyrz[M1].astype(float).mean(axis=1)
gyrz['gyrz_pMHC SEM'] = gyrz[M1].astype(float).sem(axis=1)

gyrz['ndx'] = gyrz.index

mx_gyrz_pMHC = gyrz['gyrz_pMHC'].idxmax().right

#SASA
t, t1, t2, TCR, TCR1, TCR2 = [],[],[],[],[],[]
##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            TCR.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            TCR1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            TCR2.append(float(cols[1]))


#Bin SASA  
cut_bins=np.linspace(100,300,2000).tolist()

sasa={'TCR (1ns)': TCR}
df = pd.DataFrame(sasa)
s=pd.cut(df['TCR (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'TCR (2ns)': TCR1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['TCR (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'TCR (3ns)': TCR2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['TCR (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['TCR (2ns)'] = s1['TCR (2ns)']
s2 = s2.to_frame()
SASA['TCR (3ns)'] = s2['TCR (3ns)']

l =np.linspace(100,300,2000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['TCR (1ns)','TCR (2ns)','TCR (3ns)']
SASA['TCR'] = SASA[M1].astype(float).mean(axis=1)
SASA['TCR SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_TCR = SASA['TCR'].idxmax().right
#SASA
t, t1, t2, pMHC, pMHC1, pMHC2 = [],[],[],[],[],[]
##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/SASA_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            pMHC.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/SASA_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            pMHC1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/SASA_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            pMHC2.append(float(cols[1]))

#Bin SASA  
cut_bins=np.linspace(100,300,2000).tolist()

sasa={'pMHC (1ns)': pMHC}
df = pd.DataFrame(sasa)
s=pd.cut(df['pMHC (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa1={'pMHC (2ns)': pMHC1}
df1 = pd.DataFrame(sasa1)
s1=pd.cut(df1['pMHC (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

sasa2={'pMHC (3ns)': pMHC2}
df2 = pd.DataFrame(sasa2)
s2=pd.cut(df2['pMHC (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

SASA = s.to_frame()
s1 = s1.to_frame()
SASA['pMHC (2ns)'] = s1['pMHC (2ns)']
s2 = s2.to_frame()
SASA['pMHC (3ns)'] = s2['pMHC (3ns)']

l =np.linspace(100,300,2000).tolist()
x = pd.Series(l[:-1])
SASA['SASA']= x.values

M1 = ['pMHC (1ns)','pMHC (2ns)','pMHC (3ns)']
SASA['pMHC'] = SASA[M1].astype(float).mean(axis=1)
SASA['pMHC SEM'] = SASA[M1].astype(float).sem(axis=1)

SASA['ndx'] = SASA.index

mx_pMHC = SASA['pMHC'].idxmax().right

#rmsf
t, t1, t2, TCR, TCR1, TCR2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            TCR.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            TCR1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_TCR.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            TCR2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'TCR (1ns)': TCR}
df = pd.DataFrame(rmsf)
s=pd.cut(df['TCR (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'TCR (2ns)': TCR1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['TCR (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'TCR (3ns)': TCR2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['TCR (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['TCR (2ns)'] = s1['TCR (2ns)']
s2 = s2.to_frame()
rmsf['TCR (3ns)'] = s2['TCR (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['TCR (1ns)','TCR (2ns)','TCR (3ns)']
rmsf['TCR'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['TCR SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_TCR = rmsf['TCR'].idxmax().right

#rmsf
t, t1, t2, pMHC, pMHC1, pMHC2 = [],[],[],[],[],[]

##Data Acquisition
# 1ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/70/rmsf_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t.append(float(cols[0]))
            pMHC.append(float(cols[1]))
# 2ns
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/60/rmsf_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t1.append(float(cols[0]))
            pMHC1.append(float(cols[1]))
# 3ns 
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/50/rmsf_pMHC.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t2.append(float(cols[0]))
            pMHC2.append(float(cols[1]))

#Bin rmsf  
cut_bins=np.linspace(0.01,0.5,100).tolist()
            
rmsf={'pMHC (1ns)': pMHC}
df = pd.DataFrame(rmsf)
s=pd.cut(df['pMHC (1ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf1={'pMHC (2ns)': pMHC1}
df1 = pd.DataFrame(rmsf1)
s1=pd.cut(df1['pMHC (2ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf2={'pMHC (3ns)': pMHC2}
df2 = pd.DataFrame(rmsf2)
s2=pd.cut(df2['pMHC (3ns)'], bins=cut_bins).value_counts(sort=False, normalize=True)

rmsf = s.to_frame()
s1 = s1.to_frame()
rmsf['pMHC (2ns)'] = s1['pMHC (2ns)']
s2 = s2.to_frame()
rmsf['pMHC (3ns)'] = s2['pMHC (3ns)']

l =np.linspace(0.01,0.5,100).tolist()
x = pd.Series(l[:-1])
rmsf['rmsf']= x.values

M1 = ['pMHC (1ns)','pMHC (2ns)','pMHC (3ns)']
rmsf['pMHC'] = rmsf[M1].astype(float).mean(axis=1)
rmsf['pMHC SEM'] = rmsf[M1].astype(float).sem(axis=1)

rmsf['ndx'] = rmsf.index

mx_rmsf_pMHC = rmsf['pMHC'].idxmax().right



#Compile Data for RNN

datar_labels = {'8S':resi}
datar = pd.DataFrame(datar_labels)
print(datar)
datar = datar.replace(['A','R','N','D','C','Q','E','G','H','I','L',
                         'K','M','F','P','S','T','W','Y','V']
                       ,['1','2','3','4','5','6','7','8','9','10',
                       '11','12','13','14','15','16','17','18','19','20'])

datar = datar.T


datar['Total H-Bonds'] = hb
datar['Total Contacts'] = contacts
datar['Instant H-Bonds'] = mx_insth
datar['Instant Contacts'] = mx_instc
datar['Bond Lifetime'] = time
datar['RXN Distance'] = zero
datar['Max Frequency Distance'] = mx_bndx
datar['RMSF_pMHC']=mx_rmsf_pMHC
datar['RMSF_TCR']=mx_rmsf_TCR
datar['SASA_pMHC']=mx_pMHC
datar['SASA_TCR']=mx_TCR
datar['GYR_pMHC']=mx_gyr_pMHC
datar['GYRx_pMHC']=mx_gyrx_pMHC
datar['GYRy_pMHC']=mx_gyry_pMHC
datar['GYRz_pMHC']=mx_gyrz_pMHC
datar['GYR_TCR']=mx_gyr_TCR
datar['GYRx_TCR']=mx_gyrx_TCR
datar['GYRy_TCR']=mx_gyry_TCR
datar['GYRz_TCR']=mx_gyrz_TCR
    
datar['CDR3 Distance'] = mx_CDR3dis
datar['SASA_pep'] = mx_pep
datar['SASA_CDR1a'] = mx_CDR1a
datar['SASA_CDR2a'] = mx_CDR2a
datar['SASA_CDR3a'] = mx_CDR3a
datar['SASA_CDR1b'] = mx_CDR1b
datar['SASA_CDR2b'] = mx_CDR2b
datar['SASA_CDR3b'] = mx_CDR3b
datar['SASA_MHCa'] = mx_MHCa
datar['SASA_MHCb'] = mx_MHCb
datar['RMSF_pep'] = mx_rmsf_pep
datar['RMSF_CDR1a'] = mx_rmsf_CDR1a
datar['RMSF_CDR2a'] = mx_rmsf_CDR2a
datar['RMSF_CDR3a'] = mx_rmsf_CDR3a
datar['RMSF_CDR1b'] = mx_rmsf_CDR1b
datar['RMSF_CDR2b'] = mx_rmsf_CDR2b
datar['RMSF_CDR3b'] = mx_rmsf_CDR3b
datar['RMSF_MHCa'] = mx_rmsf_MHCa
datar['RMSF_MHCb'] = mx_rmsf_MHCb
datar['GYR_pep'] = mx_gyr_pep
datar['GYRx_pep'] = mx_gyrx_pep
datar['GYRy_pep'] = mx_gyry_pep
datar['GYRz_pep'] = mx_gyrz_pep
datar['GYR_CDR1a'] = mx_gyr_CDR1a
datar['GYRx_CDR1a'] = mx_gyrx_CDR1a
datar['GYRy_CDR1a'] = mx_gyry_CDR1a
datar['GYRz_CDR1a'] = mx_gyrz_CDR1a
datar['GYR_CDR2a'] = mx_gyr_CDR2a
datar['GYRx_CDR2a'] = mx_gyrx_CDR2a
datar['GYRy_CDR2a'] = mx_gyry_CDR2a
datar['GYRz_CDR2a'] = mx_gyrz_CDR2a
datar['GYR_CDR3a'] = mx_gyr_CDR3a
datar['GYRx_CDR3a'] = mx_gyrx_CDR3a
datar['GYRy_CDR3a'] = mx_gyry_CDR3a
datar['GYRz_CDR3a'] = mx_gyrz_CDR3a
datar['GYR_CDR1b'] = mx_gyr_CDR1b
datar['GYRx_CDR1b'] = mx_gyrx_CDR1b
datar['GYRy_CDR1b'] = mx_gyry_CDR1b
datar['GYRz_CDR1b'] = mx_gyrz_CDR1b
datar['GYR_CDR2b'] = mx_gyr_CDR2b
datar['GYRx_CDR2b'] = mx_gyrx_CDR2b
datar['GYRy_CDR2b'] = mx_gyry_CDR2b
datar['GYRz_CDR2b'] = mx_gyrz_CDR2b
datar['GYR_CDR3b'] = mx_gyr_CDR3b
datar['GYRx_CDR3b'] = mx_gyrx_CDR3b
datar['GYRy_CDR3b'] = mx_gyry_CDR3b
datar['GYRz_CDR3b'] = mx_gyrz_CDR3b
datar['GYR_MHCa'] = mx_gyr_MHCa
datar['GYRx_MHCa'] = mx_gyrx_MHCa
datar['GYRy_MHCa'] = mx_gyry_MHCa
datar['GYRz_MHCa'] = mx_gyrz_MHCa
datar['GYR_MHCb'] = mx_gyr_MHCb
datar['GYRx_MHCb'] = mx_gyrx_MHCb
datar['GYRy_MHCb'] = mx_gyry_MHCb
datar['GYRz_MHCb'] = mx_gyrz_MHCb
datar['Instant Contacts CDR1a_MHCa']=mx_instc_CDR1a_MHCa
datar['Instant Contacts CDR1a_pep']=mx_instc_CDR1a_pep
datar['Total Contacts CDR1a_MHCa']=contacts_CDR1a_MHCa
datar['Total Contacts CDR1a_pep']=contacts_CDR1a_pep
datar['Instant Contacts CDR2a_MHCa']=mx_instc_CDR2a_MHCa
datar['Instant Contacts CDR2a_pep']=mx_instc_CDR2a_pep
datar['Total Contacts CDR2a_MHCa']=contacts_CDR2a_MHCa
datar['Total Contacts CDR2a_pep']=contacts_CDR2a_pep
datar['Instant Contacts CDR3a_pep']=mx_instc_CDR3a_pep
datar['Instant Contacts CDR3a_MHCa']=mx_instc_CDR3a_MHCa
datar['Total Contacts CDR3a_pep']=contacts_CDR3a_pep
datar['Total Contacts CDR3a_MHCa']=contacts_CDR3a_MHCa
datar['Instant Contacts CDR1b_MHCb']=mx_instc_CDR1b_MHCb
datar['Instant Contacts CDR1b_pep']=mx_instc_CDR1b_pep
datar['Total Contacts CDR1b_MHCb']=contacts_CDR1b_MHCb
datar['Total Contacts CDR1b_pep']=contacts_CDR1b_pep
datar['Instant Contacts CDR2b_MHCb']=mx_instc_CDR2b_MHCb
datar['Instant Contacts CDR2b_pep']=mx_instc_CDR2b_pep
datar['Total Contacts CDR2b_MHCb']=contacts_CDR2b_MHCb
datar['Total Contacts CDR2b_pep']=contacts_CDR2b_pep
datar['Instant Contacts CDR3b_pep']=mx_instc_CDR3b_pep
datar['Instant Contacts CDR3b_MHCb']=mx_instc_CDR3b_MHCb
datar['Total Contacts CDR3b_pep']=contacts_CDR3b_pep
datar['Total Contacts CDR3b_MHCb']=contacts_CDR3b_MHCb
datar['Instant H-Bonds CDR1a_MHCa']=mx_insth_CDR1a_MHCa
datar['Instant H-Bonds CDR1a_pep']=mx_insth_CDR1a_pep
datar['Total H-Bonds CDR1a_MHCa']=h_CDR1a_MHCa
datar['Total H-Bonds CDR1a_pep']=h_CDR1a_pep
datar['Instant H-Bonds CDR2a_MHCa']=mx_insth_CDR2a_MHCa
datar['Instant H-Bonds CDR2a_pep']=mx_insth_CDR2a_pep
datar['Total H-Bonds CDR2a_MHCa']=h_CDR2a_MHCa
datar['Total H-Bonds CDR2a_pep']=h_CDR2a_pep
datar['Instant H-Bonds CDR3a_pep']=mx_insth_CDR3a_pep
datar['Instant H-Bonds CDR3a_MHCa']=mx_insth_CDR3a_MHCa
datar['Total H-Bonds CDR3a_pep']=h_CDR3a_pep
datar['Total H-Bonds CDR3a_MHCa']=h_CDR3a_MHCa
datar['Instant H-Bonds CDR1b_MHCb']=mx_insth_CDR1b_MHCb
datar['Instant H-Bonds CDR1b_pep']=mx_insth_CDR1b_pep
datar['Total H-Bonds CDR1b_MHCb']=h_CDR1b_MHCb
datar['Total H-Bonds CDR1b_pep']=h_CDR1b_pep
datar['Instant H-Bonds CDR2b_MHCb']=mx_insth_CDR2b_MHCb
datar['Instant H-Bonds CDR2b_pep']=mx_insth_CDR2b_pep
datar['Total H-Bonds CDR2b_MHCb']=h_CDR2b_MHCb
datar['Total H-Bonds CDR2b_pep']=h_CDR2b_pep
datar['Instant H-Bonds CDR3b_pep']=mx_insth_CDR3b_pep
datar['Instant H-Bonds CDR3b_MHCb']=mx_insth_CDR3b_MHCb
datar['Total H-Bonds CDR3b_pep']=h_CDR3b_pep
datar['Total H-Bonds CDR3b_MHCb']=h_CDR3b_MHCb

print(datar)

datar.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/8S/fnc_bio.csv')

