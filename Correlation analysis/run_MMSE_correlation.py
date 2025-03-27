import numpy as np
import pandas as pd
import glob
import re
from utils import *
import scipy.stats
from datetime import datetime
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns


#new_data = pd.read_pickle('./data/ADNI_OASIS.pkl')
new_data = pd.read_csv('/autofs/space/genesis_001/users/Projects/PPMI/code/matlab_code/DM_OASISADNI/base_OASISADNI.csv')
new_data['DM_C'] = [i.split('connMat_Hough.mat')[0] + 'baseline_ADNIOASIS/DM_C.mat' for i in
                                 list(new_data['conn_mat'].values)]
new_data['DM_aC'] = [i.split('connMat_Hough.mat')[0] + 'baseline_ADNIOASIS/DM_aC.mat' for i in
                                 list(new_data['conn_mat'].values)]
C = read_sc(new_data, 'C')  # without log transformation
aC = read_sc(new_data,'aC') # without log transformation


DM_C = np.zeros((len(new_data), int(85 * 84 / 2)))
DM_aC = np.zeros((len(new_data), int(85 * 84 / 2)))
for i in range(len(new_data)):
    tmp = scipy.io.loadmat(new_data.loc[i, 'DM_aC'])['DM_aC']
    DM_aC[i, :] = tmp[np.triu_indices(85, k=1)]
    tmp2 = scipy.io.loadmat(new_data.loc[i, 'DM_C'])['DM_C']
    DM_C[i, :] = tmp2[np.triu_indices(85, k=1)]

nonZero_col = np.where(C.any(axis=0))[0]
zero_col = np.where(~C.any(axis=0))[0]
sc= np.delete(C, zero_col, axis=1)
sc2 = np.log(sc+1)

#new_data = new_data.rename(columns={'M/F': 'Sex', 'Site': 'SITE'})

# Combat harmonization with age as input, on log space
model, combat_age_C = harmonize(sc2, new_data, 0)
new_temp = np.zeros((len(new_data), 3570))
new_temp[:, nonZero_col] = combat_age_C
com_C = new_temp


# Combat harmonization without age as input, on log space
model, combat_noAge_C = harmonize_noAge(sc2, new_data)
new_temp = np.zeros((len(new_data), 3570))
new_temp[:, nonZero_col] = combat_noAge_C
com_noAge_C = new_temp

# Combat harmonization without age or sex as input, on log space
model, combat_noAge_noSex_C = harmonize_noAge_noSex(sc2, new_data)
new_temp = np.zeros((len(new_data), 3570))
new_temp[:, nonZero_col] = combat_noAge_noSex_C
com_noAge_noSex_C = new_temp


# Combat harmonization with age as input on aC, on log space
model, com_aC = harmonize(np.log(aC+1), new_data, 0)


# Combat harmonization without age as input on aC, on log space
model, com_noAge_aC = harmonize_noAge(np.log(aC+1), new_data)

model, com_noAge_noSex_aC = harmonize_noAge_noSex(np.log(aC+1), new_data)



#Covbat harmonization 

# In log(C+1) space
# SC = {'C': np.log(C+1), 'DM_C': np.log(DM_C+1), 'Com_C': com_C, 'Com_noAge_C': com_noAge_C, 'aC': np.log(aC+1),
#       'DM_aC': np.log(DM_aC+1),  'Com_aC':com_aC, 'Com_noAge_aC':com_noAge_aC}

### in the Iman's normalization space
SC = {'C': norm_fea(C), 'DM_C': norm_fea(DM_C), 'Com_C': norm_fea(np.exp(com_C)-1), 'Com_noAge_C': norm_fea(np.exp(com_noAge_C)-1), 
      'Com_noAge_noSex_C': norm_fea(np.exp(com_noAge_noSex_C)-1),'aC': norm_fea(aC), 'DM_aC': norm_fea(DM_aC),  'Com_aC':norm_fea(np.exp(com_aC)-1), 
      'Com_noAge_aC': norm_fea(np.exp(com_noAge_aC)-1),'Com_noAge_noSex_aC': norm_fea(np.exp(com_noAge_noSex_aC)-1),}


MMSE_index = new_data.loc[~new_data.MMSE.isna()].index
MMSE_new = new_data.loc[MMSE_index, 'MMSE'].values
age = new_data.Age.values

column_name = [i+'_MMSE' for i in list(SC.keys())] + [i+'_Age' for i in list(SC.keys())]

r = np.zeros((aC.shape[1], 20))
p = np.zeros((aC.shape[1], 20))

for j,c in enumerate(list(SC.keys())):
    fea = SC[c][MMSE_index, :]
    for i in range(aC.shape[1]):
        r[i, j], p[i, j] = scipy.stats.pearsonr(MMSE_new, fea[:, i])

for j,c in enumerate(list(SC.keys())):
    fea = SC[c]
    for i in range(aC.shape[1]):
        r[i, j+10], p[i, j+10] = scipy.stats.pearsonr(age, fea[:, i])

p = p*3570
df_corr = pd.DataFrame(data=p, columns=column_name)
#df_corr.to_pickle('/autofs/space/genesis_001/users/Projects/PPMI/code/Result_Figure/p_OASISADNI.pkl')

R = np.abs(r)
mean_r= np.nanmean(R, axis=0)
std_r = np.nanstd(R, axis=0)


# np.save('/autofs/space/genesis_001/users/Projects/PPMI/code/Result_Figure/r_OASISADNI.npy', r)
# np.save('/autofs/space/genesis_001/users/Projects/PPMI/code/Result_Figure/ab_r_OASISADNI.npy', R)

from scipy.stats import wilcoxon
for i in [0,5,10,15]:
    for j in range(i+1,i+5):
        res = wilcoxon(R[:,j], R[:,i], alternative='greater',nan_policy='omit')
        print(res)

print((df_corr < 0.05).sum())
print(df_corr.min(axis=0).apply(lambda x: f"{x:.0e}"))

for i in range(20):
    round_mean = "{:.2f}".format(mean_r[i])
    round_std = "{:.2f}".format(std_r[i])
    print( f'{round_mean} ({round_std})') 


##### correlation between C/aC and MMSE within ADNI-2 or OASIS-3 (before harmonization, individually)
## for OASIS-3

#OASIS_id = new_data.loc[new_data['SITE']=='OASIS-3'].index.tolist()
# c='C' # 'aC' or 'C'
# dataset = 'ADNI2' # 'ADNI2' or 'OASIS-3'

for c in ['C', 'aC']:
    for dataset in ['OASIS-3', 'ADNI2']:
        r = np.zeros(aC.shape[1])
        p = np.zeros(aC.shape[1])
        MMSE_index = new_data.loc[(~new_data.MMSE.isna()) & (new_data['SITE']==dataset)].index
        MMSE_OASIS = new_data.loc[MMSE_index,'MMSE'].values
        fea = SC[c][MMSE_index, :]
        for i in range(fea.shape[1]):
            r[i], p[i] = scipy.stats.pearsonr(MMSE_OASIS, fea[:, i])

        R = np.abs(r)
        mean_r= np.nanmean(R, axis=0)
        std_r = np.nanstd(R, axis=0)
        round_mean = "{:.2f}".format(mean_r)
        round_std = "{:.2f}".format(std_r)
        print('%s with %s in MMSE'%(c, dataset))
        print( f'{round_mean} ({round_std})') 


        p = p*3570
        print("{:.0e}".format(np.nanmin(p)))
        print((p < 0.05).sum())


##### correlation between C/aC and age within ADNI-2 or OASIS-3 (before harmonization, individually)
## for OASIS-3

#OASIS_id = new_data.loc[new_data['SITE']=='OASIS-3'].index.tolist()
# c='C' # 'aC' or 'C'
# dataset = 'OASIS-3' # 'ADNI2' or 'OASIS-3'

for c in ['C', 'aC']:
    for dataset in ['OASIS-3', 'ADNI2']:
        r = np.zeros(aC.shape[1])
        p = np.zeros(aC.shape[1])
        age_index = new_data.loc[(~new_data.Age.isna()) & (new_data['SITE']==dataset)].index
        age = new_data.loc[age_index,'Age'].values
        fea = SC[c][age_index, :]
        for i in range(fea.shape[1]):
            r[i], p[i] = scipy.stats.pearsonr(age, fea[:, i])

        R = np.abs(r)
        mean_r= np.nanmean(R, axis=0)
        std_r = np.nanstd(R, axis=0)
        round_mean = "{:.2f}".format(mean_r)
        round_std = "{:.2f}".format(std_r)
        print('%s with %s in age'%(c, dataset))
        print( f'{round_mean} ({round_std})') 


        p = p*3570
        print("{:.0e}".format(np.nanmin(p)))
        print((p < 0.05).sum())
