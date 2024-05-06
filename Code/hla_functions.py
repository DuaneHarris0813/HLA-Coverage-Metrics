#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:11:33 2024

@author: duaneharris
"""


# IMPORT PACKAGES:
from numpy import divide
from numpy import dot
from numpy import exp
from numpy import matmul
from numpy import multiply
from numpy import unique
from numpy import vectorize
from statistics import mean


# DICTIONARY FUNCTION:
def dict_value(d,k):
    return d[k]


# POSITION WEIGHTS:
def compute_position_weights():
    # ORIGINAL POSITION WEIGHTS
    pw = [0,    #1
          0,    #2
          0.1,  #3
          0.31, #4
          0.30, #5
          0.29, #6
          0.26, #7
          0.18, #8
          0]    #9
    return(divide(pw,sum(pw))) # Normalize Position Weights


# ENRICHMENT SCORES:
def compute_enrichment_scores():
    # ORIGINAL LOG ENRICHMENT SCORES:
    es = [0.127, #A
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
    
    # NORMALIZE ENRICHMENT SCORES:
    es = exp(es)
    es = divide(es,sum(es))
    
    # CREATE ENRICHMENT SCORE DICTIONARY:
    aa = list('ACDEFGHIKLMNPQRSTVWY')
    return(dict(zip(aa,es)))


# ALLELE FREQUENCY DICTIONARY:
def allele_frequency_dict(w):
    af = {}
    for k in w.iter_rows(min_row=3,min_col=12,max_col=13,values_only=True):
        if isinstance(k[0],str): af['HLA-'+k[0]] = k[1]
        else: break
    S = sum(af.values())    
    for k in af:
        af[k] = af[k]/S         
    return(af,S)
    
    
# IMMUNOGENICITY:
def compute_g(p,es,pw):
    g0 = [list(i) for i in p]
    g1 = vectorize(dict_value)(es,g0)
    return matmul(g1,pw)


# SIGMA:
def compute_sigma(b,p,es,pw):
    N = len(p)
    g = compute_g(p,es,pw)
    return dot(b,g)/N


# COVERAGE SCORES:
def compute_C(a,b,p,af,es,pw):
       
    # ALLELE FREQUENCIES:
    an = unique(a)
    fn = vectorize(dict_value)(af,a)
    fd = vectorize(dict_value)(af,an)
        
    # NUMERATOR:
    g  = compute_g(p,es,pw)       # Immunogenicity
    s  = multiply(b,g)            # Sigma
    c  = dot(fn,s)/len(unique(p)) # Coverage
    
    # RESULT:
    return c/sum(fd)
    

# INDIVIDUAL FREQUENCIES:
def compute_rho(f):
    M = len(f)
    rho = [[] for i in range(M)]
    for i in range(M):
        for j in range(M):
            if f[i] == -1:
                rho[i].append(-1)
            elif i==j:
                rho[i].append(f[i]**2)
            elif i>j:
                rho[i].append(2*f[i]*f[j])
            else:
                rho[i].append(-1)
    return rho


# INDIVIDUAL COVERAGE SCORES:
def compute_I(s):
    M = len(s)
    I = [[] for i in range(M)]
    for i in range(M):
        for j in range(M):
            if s[i] == 0:
                I[i].append(0)
            elif i==j:
                I[i].append(s[i])
            elif i>j:
                I[i].append(mean([s[i],s[j]]))
            else:
                I[i].append(0)
    return I
                
                
# COVERAGE RATIO:
def compute_F(a,p,b,af,es,pw,idp):
    MN = len(a)
    sn = 0; sd = 0
    g  = compute_g(p,es,pw)    
    for i in range(MN):
        t = af[a[i]] * b[i] * g[i]
        sd = sd + t
        if p[i] in idp:
            sn = sn + t
    F = sn/sd
    return F