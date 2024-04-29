#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 11:08:10 2024

@author: duaneharris
"""


# IMPORT PACKAGES:
from numpy import exp
from openpyxl import load_workbook
import hla_functions as hf


# DATA FILES:
fd  = '../Data/'
fdr = '../Results/'
fp  = 'Protein sequences.xlsx'
fr  = 'Gj.xlsx'


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


# NORMALIZATION OPTIONS:
en = 2 # Enrichment Normalization: 1=Rescale (Old Way) 2=Divide by Sum (New Way)


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


#%% LOAD PROTEIN SEQUENCE DATA ################################################


d  = load_workbook(fd+fp)
ps = []
for i,j in enumerate(d['Phase 1'].iter_rows(min_row=3, max_row=11, min_col=3, max_col=3, values_only=True)):
    ps.append(j[0])
    ps[i] = ps[i].replace('\n','')
del ps[3]; del ps[1]


# RETRIEVE FULL PROTEIN SEQUENCES:
nd = len(ps)
ep = [[] for i in range(nd)]
for i in range(nd):
    for j in range(len(ps[i])-8):
        ep[i].append(ps[i][j:j+9])
        
    
#%% MAIN ######################################################################    


gj  = []
gjo = [[] for i in range(nd)]
epo = [[] for i in range(nd)]
for i in range(nd):   
    gj.append(hf.compute_g(ep[i],es,pi))
    gjo[i],epo[i] = zip(*sorted(zip(gj[i],ep[i]),reverse=True))


#%% WRITE DATA TO FILE ########################################################


r = load_workbook(fdr+fr)
for i in range(nd):
    for j in range(len(gjo[i])):  
        r['Gj'].cell(4+j,2+2*i).value = epo[i][j]
        r['Gj'].cell(4+j,3+2*i).value = gjo[i][j]       
r.save(fdr+fr)