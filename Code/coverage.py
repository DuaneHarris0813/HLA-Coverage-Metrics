#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 16:18:13 2023

@author: duaneharris
"""

# THIS IS THE BRANCH FILE

# IMPORT PACKAGES:
from numpy import exp
from numpy import unique
import openpyxl as op
import hla_functions as hf


# DATA FOLDER:
fd  = '../Data/'
fdr = '../Results/'
fp  = 'Protein sequences.xlsx'
fr1 = 'Coverage Results.xlsx'
fr2 = 'Individual Results'


# DATA FILES:


fa =  'HLA-A Allele Frequencies.xlsx'               # File for Allele Frequency
fs = ['HLA-A_Ebola_Zaire_GP1_Weighted_Results.xlsx',
      'HLA-A_Ebola_Sudan_GP1_Weighted_Results.xlsx',      
      'HLA-A_Ebola_Zaire_NP_Weighted_Results.xlsx',
      'HLA-A_Ebola_Sudan_NP_Weighted_Results.xlsx',
      'HLA-A_SARS_Wuhan-Hu-1_Weighted_Results.xlsx',      
      'HLA-A_SARS_DeltaAY4_Weighted_Results.xlsx',
      'HLA-A_SARS_OmicronBA1_Weighted_Results.xlsx',
      'HLA-A_SARS_OmicronBA2_Weighted_Results.xlsx',
      'HLA-A_SARS_OmicronBA5_Weighted_Results.xlsx',
      'HLA-A_Burkholderia_HCP1_Weighted_Results.xlsx']
lb = ['Ebola GP1 (Zaire)',
      'Ebola GP1 (Sudan)',
      'Ebola NP (Zaire)',
      'Ebola NP (Sudan)',
      'SARS-CoV-2 Wuhan-Hu-1',
      'SARS-CoV-2 Delta AY.4',
      'SARS-CoV-2 Omicron BA.1',
      'SARS-CoV-2 Omicron BA.2',
      'SARS-CoV-2 Omicron BA.5', 
      'Burkholderia HCP1']


# fa =  'HLA-B Allele Frequencies.xlsx'               # File for Allele Frequency
# fs = ['HLA-B_Ebola_Zaire_GP1_Weighted_Results.xlsx',
#       'HLA-B_Ebola_Sudan_GP1_Weighted_Results.xlsx',      
#       'HLA-B_Ebola_Zaire_NP_Weighted_Results.xlsx',
#       'HLA-B_Ebola_Sudan_NP_Weighted_Results.xlsx',
#       'HLA-B_SARS_Wuhan-Hu-1_Weighted_Results.xlsx',      
#       'HLA-B_SARS_DeltaAY4_Weighted_Results.xlsx',
#       'HLA-B_SARS_OmicronBA1_Weighted_Results.xlsx',
#       'HLA-B_SARS_OmicronBA2_Weighted_Results.xlsx',
#       'HLA-B_SARS_OmicronBA5_Weighted_Results.xlsx',
#       'HLA-B_Burkholderia_HCP1_Weighted_Results.xlsx']
# lb = ['Ebola GP1 (Zaire)',
#       'Ebola GP1 (Sudan)',
#       'Ebola NP (Zaire)',
#       'Ebola NP (Sudan)',
#       'SARS-CoV-2 Wuhan-Hu-1',
#       'SARS-CoV-2 Delta AY.4',
#       'SARS-CoV-2 Omicron BA.1',
#       'SARS-CoV-2 Omicron BA.2',
#       'SARS-CoV-2 Omicron BA.5', 
#       'Burkholderia HCP1']


# fa =  'HLA-C Allele Frequencies.xlsx'               # File for Allele Frequency
# fs = ['HLA-C_Ebola_Zaire_GP1_Weighted_Results.xlsx',
#       'HLA-C_Ebola_Sudan_GP1_Weighted_Results.xlsx',      
#       'HLA-C_Ebola_Zaire_NP_Weighted_Results.xlsx',
#       'HLA-C_Ebola_Sudan_NP_Weighted_Results.xlsx',
#       'HLA-C_SARS_Wuhan-Hu-1_Weighted_Results.xlsx',      
#       'HLA-C_SARS_DeltaAY4_Weighted_Results.xlsx',
#       'HLA-C_SARS_OmicronBA1_Weighted_Results.xlsx',
#       'HLA-C_SARS_OmicronBA2_Weighted_Results.xlsx',
#       'HLA-C_SARS_OmicronBA5_Weighted_Results.xlsx',
#       'HLA-C_Burkholderia_HCP1_Weighted_Results.xlsx']
# lb = ['Ebola GP1 (Zaire)',
#       'Ebola GP1 (Sudan)',
#       'Ebola NP (Zaire)',
#       'Ebola NP (Sudan)',
#       'SARS-CoV-2 Wuhan-Hu-1',
#       'SARS-CoV-2 Delta AY.4',
#       'SARS-CoV-2 Omicron BA.1',
#       'SARS-CoV-2 Omicron BA.2',
#       'SARS-CoV-2 Omicron BA.5', 
#       'Burkholderia HCP1']


# OPTIONS:
ta  = 25 # Top Alleles (Max 25)


# NORMALIZATION OPTIONS:
en = 2 # Enrichment Normalization: 1=Rescale (Old Way) 2=Divide by Sum (New Way)
ft = 2 # Frequency Type: 1=Sum, 2=Average


#%% PARAMETERS ################################################################
    

# Log Enrichment Scores
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
    
    
# Position Importance
pi = [0,    #1
      0,    #2
      0.1,  #3
      0.31, #4
      0.3,  #5
      0.29, #6
      0.26, #7
      0.18, #8
      0]    #9


#%% NORMALIZATION #############################################################


# NORMALIZE POSITION IMPORTANCE:
S = sum(pi)
for i in range(len(pi)):
    pi[i] = pi[i]/S                 # Normalize so that sum(p[i])=1


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


a  = op.load_workbook(fd+fa)
rs = a.sheetnames
nr = len(rs)
af = [{} for i in range(nr)]
Z  = []
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
    Z.append(S)
    for k in af[i]:
        af[i][k] = af[i][k]/S
        
        
#%% LOAD PROTEIN SEQUENCE DATA ################################################


dp = op.load_workbook(fd+fp)
ps = []
for i,j in enumerate(dp['Phase 1'].iter_rows(min_row=2, max_row=11, min_col=3, max_col=3, values_only=True)):
    ps.append(j[0])
    ps[i] = ps[i].replace('\n','')


# DETERMINE NUMBER OF DISEASES:
if isinstance(fs,str): # If fs is a single string
    nd = 1             # Then the number of diseases is 1
else:                  # If fs is a list
    nd = len(fs)       # Then the number of diseases is length of list
    

# RETRIEVE FULL PROTEIN SEQUENCES:
Ep = [[] for i in range(nd)]
ep = [[] for i in range(nd)]
for i in range(nd):
    for j in range(len(ps[i])-8):
        Ep[i].append(ps[i][j:j+9])
        
        
# REORDER DISEASES:
ep[0] = Ep[1]
ep[1] = Ep[3]
ep[2] = Ep[2]
ep[3] = Ep[4]
ep[4] = Ep[5]
ep[5] = Ep[6]
ep[6] = Ep[7]
ep[7] = Ep[8]
ep[8] = Ep[9]
ep[9] = Ep[0]


#%% COMPUTE SCORES ############################################################

       
# INITIALIZE VARIABLES:      
an = [[] for j in range(nr)]

pd = [{} for i in range(nd)] # Immunogenicity Dictionary
up = [[] for i in range(nd)] # Names of Unique Peptides
N  = []                      # Number of Unique Peptides
M  = []                      # Total Allele Frequencies

f1 = [[] for j in range(nr)] # Allele Frequencies
f2 = []     

sg = [[[] for j in range(nr)] for i in range(nd)] # Sigma Values
st = [[[] for j in range(nr)] for i in range(nd)] # Sigma Values

c1 = [[] for i in range(nd)]                                           # Regional Coverage Scores
c2 = [[] for i in range(nd)]     


for i in range(nd): # Iterate Through Diseases 


    # LOAD DATA FOR DISEASE i:
    if isinstance(fs,str): d = op.load_workbook(fd+fs)   
    else: d = op.load_workbook(fd+fs[i])
    
        
    for j in range(nr): # Iterate Through Regions
    
    
        # GET DATA FOR REGION j:
        a1 = []; b1 = []; p1 = [] # Initialize Variables
        for k in d[rs[j]].iter_rows(min_row=2, values_only=True):
            if isinstance(k[0],str): # If the entry is not blank, Then
                a1.append(k[0]) # Allele Name
                b1.append(k[6]) # Binding Affinity
                p1.append(k[5]) # Peptide Name
        MN = len(a1) # Total Number of Alleles/Peptides
        
        
        # CHECK DATA FOR ERRORS:
        a2 = list(af[j].keys())[0:ta] # Allele Names for Top ta Alleles
        au = list(unique(a1))
        for k in range(len(au)):
            if au[k] not in a2:
                print('Allele Mismatch '+lb[i]+', '+rs[j]+', '+au[k])
        for k in range(len(a2)):
            if a2[k] not in au:
                print('Allele Mismatch '+lb[i]+', '+rs[j]+', '+a2[k])  
        pu = list(unique(p1))
        for k in range(len(pu)):
            if pu[k] not in ep[i]:
                print('Epitope Mismatch, '+lb[i]+', '+rs[j]+', '+pu[k])
        for k in range(len(ep[i])):
            if ep[i][k] not in pu:
                print('Epitope Mismatch, '+lb[i]+', '+rs[j]+', '+ep[i][k])
                
                
        # HLA-C Australia has only 22 alleles. If ta>22, then we will need to 
        # append blank entries onto a2 so that the rest of the code will run.
        if fa[4]=='C' and j==0 and len(a2)<ta: 
            for k in range(len(a2),ta):
                a2.append('')
        an[j] = a2
        
        
        # COMPUTE ALLELE AND INDIVIDUAL FREQUENCIES:
        if i==0:
            for k in range(ta):
                if a2[k]=='':
                    f1[j].append(-1)
                else:
                    f1[j].append(af[j][a2[k]])
            f2.append(hf.compute_rho(f1[j]))
            for k in range(ta):
                if f1[j][k]==-1: f1[j][k]=0
                
                
        # COMPUTE SIGMA FOR EACH ALLELE:
        for k in range(ta):
            if a2[k] == '':
                sg[i][j].append(0)
            else:
                ia = [l for l in range(MN) if a1[l] == a2[k]]
                sg[i][j].append(hf.compute_sigma([b1[l] for l in ia],[p1[l] for l in ia],es,pi))
                
                
        # COMPUTE REGIONAL COVERAGE SCORE:
        c1[i].append(hf.compute_C(a1,b1,p1,af[j],es,pi))
        
        
        # COMPUTE INDIVIDUAL COVERAGE SCORES:
        c2[i].append(hf.compute_I(sg[i][j]))


#%% WRITE RESULTS TO FILE #####################################################


# LOAD RESULTS SHEETS:
r = op.load_workbook(fdr+fr1)


# ALLELE TYPE:
at = fa[4]


# COVERAGE RESULTS:
w = r['Coverage']
if at == 'A':
    rrw = 4
elif at == 'B':
    rrw = 18
elif at == 'C':
    rrw = 32
for i in range(nd):
    for j in range(nr):
        w.cell(rrw+j,i+3).value = c1[i][j]
        
        
# FREQ VS SIGMA RESULTS:
if at == 'A':
    rrw = 4
elif at == 'B':
    rrw = 33
elif at == 'C':
    rrw = 62
for j in range(nr):
    w = r[rs[j]]
    for i in range(nd):
        for k in range(ta):
            if i==0:
                w.cell(rrw+k,2).value = an[j][k]
                w.cell(rrw+k,3).value = f1[j][k]
            w.cell(rrw+k,4+i).value = sg[i][j][k]
        
        
# SAVE RESULTS 1:
r.save(fdr+fr1)


# INDIVIDUAL RESULTS:
for j in range(nr):   
    
    # LOAD RESULT FILE:
    r = op.load_workbook(fdr+fr2+' '+rs[j]+'.xlsx')
     
    # FREQUENCY:
    w = r['Frequency']
    for k in range(ta):
        w.cell(rrw+k,2).value = an[j][k]
        w.cell(rrw-1,3+k).value = an[j][k]
    for k in range(ta):
        for l in range(ta):
            w.cell(rrw+k,3+l).value = f2[j][k][l]
            
    # COVERAGE:
    for i in range(nd):
        w = r[lb[i]]
        for k in range(ta):
            w.cell(rrw+k,2).value = an[j][k]
            w.cell(rrw-1,3+k).value = an[j][k]
        for k in range(ta):
            for l in range(ta):
                w.cell(rrw+k,3+l).value = c2[i][j][k][l]
    
    # SAVE RESULTS 2
    r.save(fdr+fr2+' '+rs[j]+'.xlsx')        
