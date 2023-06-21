#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import lightgbm as lgb
import glob
import sys

thr = 0.6112519174331289


#with open(sys.argv[1], "r") as file:
#    labels = eval(file.readline())


data = pd.read_csv(sys.argv[1], index_col=0)
data

models = glob.glob('MODELS/*.txt')
main_columns = ['GAtotal', 'BWg', 'SGA3rd', 'TempC', 'BEmmolL', 'CRIBII', 'Sex_F',
       'Race_Black', 'Race_White', 'VentDaySIMV_CPAP', 'FiO2Day3', 'FiO27',
       'Surfactantdoses', 'CRPpeak', 'CRP2', 'WBC4', 'WBC', 'Bands10', 'Bands',
       'Metamyeloctye', 'Neutrophils', 'AbxDays', 'abxdoses', 'IT', 'IT02',
       'BMIkgm2', 'MaternalAge', 'Primagravida', 'AntenatalSteroids',
       'Prenatalabx72h', 'VentDaySIMV_Cannula', 'VentDaySIMV_HFV',
       'VentDaySIMV_SIMV', 'VentDay7_CPAP', 'VentDay7_Cannula', 'VentDay7_HFV',
       'VentDay7_SIMV', 'Race_Asian', 'Bacteria;Actinobacteriota',
       'Bacteria;Proteobacteria', 'Bacteria;Bacteroidota',
       'Bacteria;Fusobacteriota', 'Bacteria;Verrucomicrobiota',
       'Bacteria;Patescibacteria', 'Bacteria;Campylobacterota',
       'Bacteria;Bdellovibrionota', 'Bacteria;Desulfobacterota',
       'Bacteria;Spirochaetota', 'Bacteria;Synergistota', 'Fungi;Ascomycota',
       'Fungi;Basidiomycota', 'Fungi;Mortierellomycota', 'Fungi;Rozellomycota',
       'Fungi;Mucoromycota', 'Bacteria;Firmicutes']

final_vars =['BWg',
 'CRIBII',
 'CRPpeak',
 'Fungi;Mortierellomycota',
 'FiO27',
 'Fungi;Rozellomycota',
 'Bacteria;Bacteroidota',
 'Fungi;Ascomycota',
 'Fungi;Basidiomycota',
 'FiO2Day3',
 'IT',
 'GAtotal',
 'Bacteria;Proteobacteria',
 'BMIkgm2',
 'MaternalAge',
 'Primagravida']

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

