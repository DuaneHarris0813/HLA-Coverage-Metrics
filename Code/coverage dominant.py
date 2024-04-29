#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:00:19 2024

@author: duaneharris
"""

# IMPORT PACKAGES:
from numpy import exp
from openpyxl import load_workbook
import hla_functions as hf


# DATA FOLDER:
fd  = '../Data/'
fdr = '../Results/'
fr  = 'Coverage Results.xlsx'


# DATA FILES:
fa =  'HLA-A Allele Frequencies.xlsx'               # File for Allele Frequency
fs = ['HLA-A_Ebola_Zaire_GP1_Weighted_Results.xlsx',
      'HLA-A_Ebola_Sudan_GP1_Weighted_Results.xlsx',      
      'HLA-A_SARS_Wuhan-Hu-1_Weighted_Results.xlsx',     
      'HLA-A_SARS_DeltaAY4_Weighted_Results.xlsx',
      'HLA-A_SARS_OmicronBA1_Weighted_Results.xlsx',
      'HLA-A_SARS_OmicronBA2_Weighted_Results.xlsx',
      'HLA-A_SARS_OmicronBA5_Weighted_Results.xlsx']
lb = ['Ebola GP1 (Zaire)',
      'Ebola GP1 (Sudan)',
      'SARS-CoV-2 Wuhan-Hu-1',
      'SARS-CoV-2 Delta AY.4',
      'SARS-CoV-2 Omicron BA.1',
      'SARS-CoV-2 Omicron BA.2',
      'SARS-CoV-2 Omicron BA.5']


# fa =  'HLA-B Allele Frequencies.xlsx'               # File for Allele Frequency
# fs = ['HLA-B_Ebola_Zaire_GP1_Weighted_Results.xlsx',
#       'HLA-B_Ebola_Sudan_GP1_Weighted_Results.xlsx',      
#       'HLA-B_SARS_Wuhan-Hu-1_Weighted_Results.xlsx',      
#       'HLA-B_SARS_DeltaAY4_Weighted_Results.xlsx',
#       'HLA-B_SARS_OmicronBA1_Weighted_Results.xlsx',
#       'HLA-B_SARS_OmicronBA2_Weighted_Results.xlsx',
#       'HLA-B_SARS_OmicronBA5_Weighted_Results.xlsx']
# lb = ['Ebola GP1 (Zaire)',
#       'Ebola GP1 (Sudan)',
#       'SARS-CoV-2 Wuhan-Hu-1',
#       'SARS-CoV-2 Delta AY.4',
#       'SARS-CoV-2 Omicron BA.1',
#       'SARS-CoV-2 Omicron BA.2',
#       'SARS-CoV-2 Omicron BA.5']


# fa =  'HLA-C Allele Frequencies.xlsx'               # File for Allele Frequency
# fs = ['HLA-C_Ebola_Zaire_GP1_Weighted_Results.xlsx',
#       'HLA-C_Ebola_Sudan_GP1_Weighted_Results.xlsx',      
#       'HLA-C_SARS_Wuhan-Hu-1_Weighted_Results.xlsx',      
#       'HLA-C_SARS_DeltaAY4_Weighted_Results.xlsx',
#       'HLA-C_SARS_OmicronBA1_Weighted_Results.xlsx',
#       'HLA-C_SARS_OmicronBA2_Weighted_Results.xlsx',
#       'HLA-C_SARS_OmicronBA5_Weighted_Results.xlsx']
# lb = ['Ebola GP1 (Zaire)',
#       'Ebola GP1 (Sudan)',
#       'SARS-CoV-2 Wuhan-Hu-1',
#       'SARS-CoV-2 Delta AY.4',
#       'SARS-CoV-2 Omicron BA.1',
#       'SARS-CoV-2 Omicron BA.2',
#       'SARS-CoV-2 Omicron BA.5']


# NORMALIZATION OPTIONS:
en = 2 # Enrichment Normalization: 1=Rescale (Old Way) 2=Divide by Sum (New Way)
ft = 2 # Frequency Type: 1=Sum, 2=Average


#%% PARAMETERS ################################################################
    

# REGIONS:
rg = ['Australia', 
      'Europe',
      'North Africa',
      'North America',
      'North East Asia',
      'Oceania',
      'South and Central America',
      'South Asia',
      'South East Asia',
      'Sub-Saharan Africa',
      'Western Asia']


# LOG ENRICHMENT SCORES:
le = [0.127, #A
     -0.175, #C
      0.072, #D
      0.325, #E
      0.380, #F
      0.110, #G
      0.105, #H
      0.432, #I
     -0.700, #K
     -0.036, #L
     -0.570, #M
     -0.021, #N
     -0.036, #P
     -0.376, #Q
      0.168, #R
     -0.537, #S
      0.126, #T
      0.134, #V
      0.719, #W
     -0.012] #Y
    
    
# POSITION IMPORTANCE
pw = [0,    #1
      0,    #2
      0.1,  #3
      0.31, #4
      0.3,  #5
      0.29, #6
      0.26, #7
      0.18, #8
      0]    #9


idpe = ['ATDVPSATK',
        'TDVPSATKR',
        'GFRSGVPPK',
        'AENCYNLEI',
        'RLASTVIYR',
        'TEDPSSGYY',
        'DTTIGEWAF',
        'TTIGEWAFW',
        'NQDGLICGL',
        'TELRTFSIL',
        'ALFCICKFV',
        'LFCICKFVF']
idps = ['YLQPRTFLL',
        'GVYFASTEK',
        'NLNESLIDL',
        'TLDSKTQSL',
        'RLQSLQTYV',
        'QIYKTPPIK']


#%% NORMALIZATION #############################################################


# NORMALIZE POSITION IMPORTANCE:
S = sum(pw)
for i in range(len(pw)):
    pw[i] = pw[i]/S                 # Normalize so that sum(p[i])=1


# NORMALIZE ENRICHMENT SCORES:
for i in range(len(le)):
    le[i] = exp(le[i])              # Undo Natural Log
M = max(le); m = min(le)            # M = Max Enrichment, m = Min Enrichment
S = sum(le)
for i in range(len(le)):        
    if en==1:       
        le[i] = (le[i] - m) / (M-m) # Rescale Enrichment Scores
    if en==2:
        le[i] = le[i]/S             # Normalize Enrichment Scores
        
        
# CREATE ENRICHMENT SCORE DICTIONARY:
aa = list('ACDEFGHIKLMNPQRSTVWY')
es = dict(zip(aa,le))


#%% LOAD ALLELE FREQUENCY DATA ################################################


a  = load_workbook(fd+fa)
rs = a.sheetnames
nr = len(rs)
af = [{} for i in range(nr)]
for i in range(nr):
    ad = []
    if ft==1: m1=18; m2=19
    if ft==2: m1=12; m2=13
    for k in a[rs[i]].iter_rows(min_row=3, min_col=m1, max_col=m2, values_only=True):
        if isinstance(k[0], str): ad.append(k)
        else: break
    na = len(ad)
    S = 0
    for k in range(na):
        af[i]['HLA-'+ad[k][0]]=ad[k][1]
        S = S + ad[k][1]
    for k in af[i]:
        af[i][k] = af[i][k]/S
        
        
#%% COMPUTE SCORES ############################################################


# DETERMINE NUMBER OF DISEASES:
if isinstance(fs,str): # If fs is a single string
    nd = 1             # Then the number of diseases is 1
else:                  # If fs is a list
    nd = len(fs)       # Then the number of diseases is length of list

       
# INITIALIZE COVERAGE VARIABLE:      
F = [[] for i in range(nd)]
 

# COMPUTE COVERAGE RATIOS:
for i in range(nd): # Iterate Through Diseases 


    # LOAD DATA FOR DISEASE i:
    if isinstance(fs,str): d = load_workbook(fd+fs)   
    else: d = load_workbook(fd+fs[i])
    
    
    # READ NECESSARY DATA:
    if i<=1: idp=idpe
    else: idp=idps
    
        
    for j in range(nr): # Iterate Through Regions
    
    
        # GET DATA FOR REGION j:
        a1 = []; p1 = []; b1 = [];
        for k in d[rg[j]].iter_rows(min_row=2, min_col=1, max_col=7, values_only=True):
            a1.append(k[0])     # Allele for Disease i Region j
            p1.append(k[5])     # Peptides for Disease i Region j
            b1.append(k[6])     # Binding Affinities for Disease i Region j
   
   
        # COMPUTE F RATIO:     
        F[i].append(hf.compute_F(a1,p1,b1,af[j],es,pw,idp))
        
        
#%% SAVE RESULTS:


# LOAD RESULTS SHEET:
r = load_workbook(fdr+fr)
w = r['Dominant Coverage']


# DETERMINE ALLELE TYPE:
at = fa[4]
if at=='A':
    rrw = 4
elif at=='B':
    rrw = 18
elif at=='C':
    rrw = 32
    
    
# WRITE TO FILE:
for i in range(nd):
    for j in range(nr):
        w.cell(rrw+j,i+3).value = F[i][j]


# SAVE FILE:
r.save(fdr+fr)