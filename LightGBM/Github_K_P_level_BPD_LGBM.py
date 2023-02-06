#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import lightgbm as lgb
import glob
import sys

thr = 0.3565747661240787


#with open(sys.argv[1], "r") as file:
#    labels = eval(file.readline())


data = pd.read_csv(sys.argv[1], index_col=0)
data

models = glob.glob('Models/*.txt')

main_columns = ['GAtotal', 'BWg', 'SGA3rd', 'TempC',
                              'Sex_F', 'Race_Black', 'Race_White',
                               'VentDay1_CPAP', 'VentDay1_Hood', 'VentDay1_SIMV', 'Unassigned',
       'Archaea', 'Bacteria', 'Fungi', 'BacteriaThermi',
       'BacteriaFusobacteria', 'BacteriaAcidobacteria',
       'BacteriaActinobacteria', 'FungiAscomycota', 'BacteriaTM6',
       'ArchaeaEuryarchaeota', 'Fungiunidentified', 'ArchaeaParvarchaeota',
       'BacteriaTenericutes', 'BacteriaSpirochaetes', 'BacteriaProteobacteria',
       'BacteriaFibrobacteres', 'BacteriaDeferribacteres',
       'BacteriaVerrucomicrobia', 'FungiMucoromycota', 'BacteriaTM7',
       'FungiBasidiomycota', 'BacteriaChloroflexi', 'BacteriaCyanobacteria',
       'BacteriaPlanctomycetes', 'BacteriaBacteroidetes', 'BacteriaFirmicutes',
       'BacteriaLentisphaerae', 'BacteriaGN04', 'BacteriaArmatimonadetes',
       'BacteriaGemmatimonadetes', 'FungiMortierellomycota',
       'ArchaeaCrenarchaeota', 'BPD' ]

final_vars = ['GAtotal',
 'BWg',
 'TempC',
 'Fungi',
 'FungiAscomycota',
 'BacteriaProteobacteria',
 'BacteriaVerrucomicrobia',
 'FungiBasidiomycota']


input_frame = pd.DataFrame(0, index=np.arange(len(data)), columns=main_columns)
for var_i in final_vars:
    input_frame[var_i] = data[var_i].values



pred_prob=[]

for i in range(5):
    
    clf_name = models[i]
    clf = lgb.Booster(model_file=clf_name)
    
    pred_prob.append(clf.predict(input_frame, num_iteration=clf.best_iteration))
    


oof_hidden = np.mean(pred_prob, axis = 0) 
binary_hidden = np.array([1 if i > thr else 0 for i in oof_hidden])


with open("risk_predictions.txt", "w") as file:
    file.write(str(oof_hidden))


with open("binary_predictions.txt", "w") as file:
    file.write(str(binary_hidden))


#print('AUC is ' + str(roc_auc_score(labels, pred_prob_arg)) )

