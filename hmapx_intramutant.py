import numpy as np
import pandas as pd
import statistics as st
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import gromacs
from gromacs.formats import XPM



L1_100 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/hmap.csv', index_col = 0, header=None)
L1_95 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/130/hmap.csv', index_col = 0, header=None)
L1_90 = pd.read_csv(r'/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/120/hmap.csv', index_col = 0, header=None)

#L1
x_L1_100, x_L1_95, x_L1_90 = [],[],[]
t_L1_100, t_L1_95, t_L1_90 = [],[],[]

with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/140/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t_L1_100.append(float(cols[0]))
            x_L1_100.append(float(cols[1]))
            
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/130/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t_L1_95.append(float(cols[0]))
            x_L1_95.append(float(cols[1]))
            
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/120/pullxcf.xvg") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            t_L1_90.append(float(cols[0]))
            x_L1_90.append(float(cols[1]))


x_L1_100 = x_L1_100[::10]            
x_L1_95 = x_L1_95[::10]
x_L1_90 = x_L1_90[::10]

rows_100, cols_100 = L1_100.shape
rows_95, cols_95 = L1_95.shape
rows_90, cols_90 = L1_90.shape

cut_L1_100 = list(range(len(x_L1_100),cols_100 +1))
cut_L1_95 = list(range(len(x_L1_95), cols_95 +1))
cut_L1_90 = list(range(len(x_L1_90),cols_90 +1))

print(L1_100.shape)
print(L1_95.shape)
print(L1_90.shape)

L1_100 = L1_100.drop(columns=cut_L1_100)
L1_95 = L1_95.drop(columns=cut_L1_95)
L1_90 = L1_90.drop(columns=cut_L1_90)

print(L1_100.shape)
print(L1_95.shape)
print(L1_90.shape)

print(len(x_L1_100))
print(len(x_L1_95))
print(len(x_L1_90))

print(L1_100.shape)
print(L1_95.shape)
print(L1_90.shape)

L1_100['Total Occupancy (%)'] = L1_100.mean(axis=1).round(decimals=3)
L1_100['Total Occupancy (%)'] = 100*L1_100['Total Occupancy (%)']
L1_100_5 = L1_100[L1_100['Total Occupancy (%)']> 0]

L1_95['Total Occupancy (%)'] = L1_95.mean(axis=1).round(decimals=3)
L1_95['Total Occupancy (%)'] = 100*L1_95['Total Occupancy (%)']
L1_95_5 = L1_95[L1_95['Total Occupancy (%)']> 0]

L1_90['Total Occupancy (%)'] = L1_90.mean(axis=1).round(decimals=3)
L1_90['Total Occupancy (%)'] = 100*L1_90['Total Occupancy (%)']
L1_90_5 = L1_90[L1_90['Total Occupancy (%)']> 0]

ndx_L1100 = pd.Index(L1_100_5.index.tolist())
ndx_L195 = pd.Index(L1_95_5.index.tolist())
ndx_L190 = pd.Index(L1_90_5.index.tolist())

###L1###

#L1 Common Index
L1_100_95= ndx_L195.difference(ndx_L1100, sort=False).tolist()
L1_100_90= ndx_L190.difference(ndx_L1100, sort=False).tolist()
L1_100_diff = list(set(L1_100_95+L1_100_90))
L1_95_100 = ndx_L1100.difference(ndx_L195, sort=False).tolist()
L1_95_90 = ndx_L190.difference(ndx_L195, sort=False).tolist()
L1_95_diff = list(set(L1_95_100+L1_95_90)) 
L1_90_100 = ndx_L1100.difference(ndx_L190, sort=False).tolist()
L1_90_95 = ndx_L195.difference(ndx_L190, sort=False).tolist()
L1_90_diff = list(set(L1_90_100+L1_90_95))

##L1 Common Index Addition Loops##

#Add index differences to 100 ns
for i in range(len(L1_100_diff)):
    for j in range(len (L1_100.index.tolist())):
        if L1_100_diff[i] == L1_100.index.tolist()[j]:
            L1_100_5 = L1_100_5.append(L1_100.iloc[j])
            break
        if  j == len (L1_100.index.tolist())-1:
                add_zeros = np.zeros(shape=(1, len(L1_100.columns)-1))
                zero = pd.DataFrame(add_zeros, index = [L1_100_diff[i]])
                L1_100_5 = L1_100_5.append(zero)
L1_100_5['Total Occupancy (%)'] = 100*L1_100_5.mean(axis=1).round(decimals=3)
L1_100_5 = L1_100_5.sort_index(axis=0, ascending=True)

#Add index differences to 95 ns
for i in range(len(L1_95_diff)):
    for j in range(len (L1_95.index.tolist())):
        if L1_95_diff[i] == L1_95.index.tolist()[j]:
            L1_95_5 = L1_95_5.append(L1_95.iloc[j])
            break
        if  j == len (L1_95.index.tolist())-1:
            add_zeros = np.zeros(shape=(1, len(L1_95.columns)-1))
            zero = pd.DataFrame(add_zeros, index = [L1_95_diff[i]])
            L1_95_5 = L1_95_5.append(zero)
    
L1_95_5['Total Occupancy (%)'] = 100*L1_95_5.mean(axis=1).round(decimals=3)
L1_95_5 = L1_95_5.sort_index(axis=0, ascending=True)

#Add index differences to 90 ns
for i in range(len(L1_90_diff)):
    for j in range(len (L1_90.index.tolist())):
        if L1_90_diff[i] == L1_90.index.tolist()[j]:
            L1_90_5 = L1_90_5.append(L1_90.iloc[j])
            break
        if  j == len (L1_90.index.tolist())-1:
                add_zeros = np.zeros(shape=(1, len(L1_90.columns)-1))
                zero = pd.DataFrame(add_zeros, index = [L1_90_diff[i]])
                L1_90_5 = L1_90_5.append(zero)           
L1_90_5['Total Occupancy (%)'] = 100*L1_90_5.mean(axis=1).round(decimals=3)
L1_90_5 = L1_90_5.sort_index(axis=0, ascending=True)

###L1###

FS=7
bins=np.linspace(5.884,9.000,200).tolist()
xticks=bins[::20]
xticks=np.around(xticks,2)
                
#Convert to COM Distance & Average Mutant Maps
print(L1_100_5.shape)
print(L1_95_5.shape)
print(L1_90_5.shape)

L1_100_ocpy = pd.DataFrame()
L1_100_ocpy['Total Occupancy (%)'] = L1_100_5['Total Occupancy (%)']
L1_100f = L1_100_5.drop('Total Occupancy (%)', axis=1)
#100f_rows, 100f_cols= L1_100f.shape
L1_100f.columns = x_L1_100
L1_100f=L1_100f.sort_index(axis=1, ascending=True)
L1_100fx = pd.DataFrame()
L1_100fx = L1_100fx.append(L1_100f.groupby(pd.cut(L1_100f.columns.tolist(), bins=bins),axis=1).mean())
L1_100fx = L1_100fx.fillna(0)

MHC_100 = L1_100_5[L1_100_5.index.str.contains('MHC')]
MHC_100_ocpy = pd.DataFrame()
MHC_100_ocpy['Total Occupancy (%)'] = MHC_100['Total Occupancy (%)']
MHC_100 = MHC_100.drop('Total Occupancy (%)', axis=1)
MHC_rows, MHC_cols = MHC_100.shape
#print(MHC_100)
#print(MHC_100.shape)
MHC_100.columns = x_L1_100
MHC_100 = MHC_100.sort_index(axis=1, ascending=True)
MHC_L1_100x = pd.DataFrame()
MHC_L1_100x = MHC_L1_100x.append(MHC_100.groupby(pd.cut(MHC_100.columns.tolist(), bins=bins),axis=1).mean())
MHC_L1_100x = MHC_L1_100x.fillna(0)

L1_95_ocpy = pd.DataFrame()
L1_95_ocpy['Total Occupancy (%)'] = L1_95_5['Total Occupancy (%)']
L1_95f = L1_95_5.drop('Total Occupancy (%)', axis=1)
#95f_rows, 95f_cols= L1_95f.shape
L1_95f.columns = x_L1_95
L1_95f=L1_95f.sort_index(axis=1, ascending=True)
L1_95fx = pd.DataFrame()
L1_95fx = L1_95fx.append(L1_95f.groupby(pd.cut(L1_95f.columns.tolist(), bins=bins),axis=1).mean())
L1_95fx = L1_95fx.fillna(0)
#print(L1_95fx)

MHC_95 = L1_95_5[L1_95_5.index.str.contains('MHC')]
MHC_95_ocpy = pd.DataFrame()
MHC_95_ocpy['Total Occupancy (%)'] = MHC_95['Total Occupancy (%)']
MHC_95 = MHC_95.drop('Total Occupancy (%)', axis=1)
#MHC_95 = MHC_95.drop(MHC_95.columns[-1],axis=1)
MHC_95 = MHC_95.fillna(0)
print(MHC_95)
print(len(x_L1_95))
MHC_95.columns = x_L1_95
MHC_95 = MHC_95.sort_index(axis=1, ascending=True)
MHC_L1_95x = pd.DataFrame()
MHC_L1_95x = MHC_L1_95x.append(MHC_95.groupby(pd.cut(MHC_95.columns.tolist(), bins=bins),axis=1).mean())
MHC_L1_95x = MHC_L1_95x.fillna(0)
print(MHC_L1_95x)

L1_90_ocpy = pd.DataFrame()
L1_90_ocpy['Total Occupancy (%)'] = L1_90_5['Total Occupancy (%)']
L1_90f = L1_90_5.drop('Total Occupancy (%)', axis=1)
#90f_rows, 90f_cols= L1_90f.shape
L1_90f.columns = x_L1_90
L1_90f=L1_90f.sort_index(axis=1, ascending=True)
L1_90fx = pd.DataFrame()
L1_90fx = L1_90fx.append(L1_90f.groupby(pd.cut(L1_90f.columns.tolist(), bins=bins),axis=1).mean())
L1_90fx = L1_90fx.fillna(0)

MHC_90 = L1_90_5[L1_90_5.index.str.contains('MHC')]
MHC_90_ocpy = pd.DataFrame()
MHC_90_ocpy['Total Occupancy (%)'] = MHC_90['Total Occupancy (%)']
MHC_90 = MHC_90.drop('Total Occupancy (%)', axis=1)
#MHC_90 = MHC_90.drop(MHC_90.columns[-1],axis=1)
MHC_90 = MHC_90.fillna(0)
MHC_90.columns = x_L1_90
MHC_90 = MHC_90.sort_index(axis=1, ascending=True)
MHC_L1_90x = pd.DataFrame()
MHC_L1_90x = MHC_L1_90x.append(MHC_90.groupby(pd.cut(MHC_90.columns.tolist(), bins=bins),axis=1).mean())
MHC_L1_90x = MHC_L1_90x.fillna(0)

MHC_AVG = pd.concat([MHC_L1_100x,MHC_L1_95x,MHC_L1_90x]).groupby(level=0).mean()
pMHC_AVG = pd.concat([L1_100fx,L1_95fx,L1_90fx]).groupby(level=0).mean()

MHC_AVG.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/hmapx_mhc.csv', header=None)
pMHC_AVG.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/bio/hCD9/hmapx_pmhc.csv', header=None)
                
##plot=ax[0].imshow(MHC_AVG,aspect='auto', cmap=plt.cm.hot, interpolation='nearest')
#ax[0].set_yticks(range(len(MHC_AVG.index)))
#ax[0].set_yticklabels(MHC_AVG.index, fontsize=FS)
#ax[0].set_title('High Occupancy Contact Heat Map (L1, AVG)',fontsize=20)


#fig.colorbar(plot, ax=ax[0:], location='right')
#fig.text(-0.03,0.5, 'Bond Pair (Donor:Acceptor)', ha="center", va="center", rotation=90, fontsize=15)
#fig.text(0.85,0.5, 'Fractional Occupancy', ha="center", va="center", rotation=270,fontsize=15)
#plt.savefig('/Users/zrollins/Box/DMF5_MART1/L1+/100/contact_occupancy_avgx.png', bbox_inches='tight', dpi=300)





