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

models = glob.glob('Models/*.txt')
main_columns = ['SampleCode', 'GAtotal', 'BWg', 'SGA3rd', 'TempC', 'BEmmolL', 'CRIBII',
       'Sex_F', 'Race_Black', 'Race_White', 'VentDaySIMV_CPAP',
       'VentDaySIMV_Cannula', 'VentDaySIMV_HFV', 'VentDaySIMV_SIMV',
       'VentDay7_CPAP', 'VentDay7_Cannula', 'VentDay7_HFV', 'VentDay7_SIMV',
       'Race_Asian', 'Unassigned;', 'k__Archaea;', 'k__Fungi;', 'k__Bacteria;',
       'k__Fungi; p__Basidiomycota;', 'k__Bacteria; p__Armatimonadetes;',
       'k__Fungi; p__Aphelidiomycota;', 'k__Bacteria; p__Proteobacteria;',
       'k__Bacteria; p__Deferribacteres;', 'k__Fungi; p__GS19;',
       'k__Bacteria; p__Chloroflexi;', 'k__Bacteria; p__WS2;',
       'k__Fungi; p__Blastocladiomycota;', 'k__Bacteria; p__SBR1093;',
       'k__Bacteria; p__Acidobacteria;', 'k__Bacteria; p__AC1;',
       'k__Bacteria; p__SC4;', 'k__Bacteria; p__Caldithrix;',
       'k__Bacteria; p__TM7;', 'k__Fungi; p__GS01;',
       'k__Fungi; p__Monoblepharomycota;', 'k__Fungi; p__unidentified;',
       'k__Bacteria; p__Planctomycetes;', 'k__Bacteria; p__;',
       'k__Bacteria; p__WWE1;', 'k__Bacteria; p__Synergistetes;',
       'k__Bacteria; p__Tenericutes;', 'k__Bacteria; p__Firmicutes;',
       'k__Bacteria; p__GN04;', 'k__Bacteria; p__Bacteroidetes;',
       'k__Bacteria; p__Nitrospirae;', 'k__Bacteria; p__GN02;',
       'k__Bacteria; p__WS3;', 'k__Bacteria; p__OP8;', 'k__Bacteria; p__OD1;',
       'k__Bacteria; p__TM6;', 'k__Fungi; p__Ascomycota;',
       'k__Bacteria; p__NKB19;', 'k__Bacteria; p__WS4;',
       'k__Bacteria; p__FCPU426;', 'k__Bacteria; p__KSB3;',
       'k__Bacteria; p__SAR406;', 'k__Fungi; p__Glomeromycota;',
       'k__Bacteria; p__Cyanobacteria;', 'k__Bacteria; p__Elusimicrobia;',
       'k__Bacteria; p__GOUTA4;', 'k__Bacteria; p__BRC1;',
       'k__Bacteria; p__Spirochaetes;', 'k__Bacteria; p__TPD-58;',
       'k__Bacteria; p__NC10;', 'k__Bacteria; p__Verrucomicrobia;',
       'k__Archaea; p__Euryarchaeota;', 'k__Bacteria; p__1Caldithrix2;',
       'k__Fungi; p__Chytridiomycota;', 'k__Bacteria; p__Gemmatimonadetes;',
       'k__Bacteria; p__SR1;', 'k__Fungi; p__Entomophthoromycota;',
       'k__Bacteria; p__Fibrobacteres;', 'k__Fungi; p__Mortierellomycota;',
       'k__Bacteria; p__WS6;', 'k__Bacteria; p__Fusobacteria;',
       'k__Bacteria; p__Chlorobi;', 'k__Bacteria; p__1Thermi2;',
       'k__Bacteria; p__Actinobacteria;', 'k__Fungi; p__Rozellomycota;',
       'k__Fungi; p__Mucoromycota;', 'k__Bacteria; p__LCP-89;',
       'k__Bacteria; p__Hyd24-12;', 'k__Archaea; p__Crenarchaeota;',
       'k__Bacteria; p__Chlamydiae;']

final_vars = ['BWg',
 'CRIBII',
 'GAtotal',
 'VentDay7_SIMV',
 'k__Fungi; p__Chytridiomycota;',
 'k__Fungi;',
 'k__Bacteria;',
 'k__Fungi; p__Ascomycota;',
 'k__Fungi; p__Rozellomycota;',
 'BEmmolL',
 'k__Bacteria; p__;',
 'k__Bacteria; p__Proteobacteria;',
 'k__Fungi; p__Mucoromycota;',
 'k__Fungi; p__Basidiomycota;',
 'VentDaySIMV_SIMV',
 'TempC']


input_frame = pd.DataFrame(0, index=np.arange(len(data)), columns=main_columns)
for var_i in final_vars:
    input_frame[var_i] = data[var_i].values



pred_prob=[]

for i in range(5):
    
    clf_name = models[i]
    clf = lgb.Booster(model_file=clf_name)
    
    pred_prob.append(clf.predict(input_frame.iloc[:,1:], num_iteration=clf.best_iteration))
    


oof_hidden = np.mean(pred_prob, axis = 0) 
binary_hidden = np.array([1 if i > thr else 0 for i in oof_hidden])


with open("risk_predictions.txt", "w") as file:
    file.write(str(oof_hidden))


with open("binary_predictions.txt", "w") as file:
    file.write(str(binary_hidden))


#print('AUC is ' + str(roc_auc_score(labels, pred_prob_arg)) )

