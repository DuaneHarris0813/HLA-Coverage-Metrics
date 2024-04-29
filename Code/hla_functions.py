#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:11:33 2024

@author: duaneharris
"""


# IMPORT PACKAGES:
from numpy import dot
from numpy import matmul
from numpy import multiply
from numpy import unique
from numpy import vectorize
from statistics import mean


def dict_value(d,k):
    return d[k]


# IMMUNOGENICITY
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


# # INDIVIDUAL COVERAGE RATIO:
# def compute_H(a,p,b,af,es,pw,idp):
    
#     # UNIQUE ALLELES:
#     an = unique(a) # Allele Names
#     M  = len(an)   # Number of Alleles
#     MN = len(a)    # Number of Alleles * Number of Epitopes
    
#     # COMPUTE INDIVIDUAL FREQUENCIES:
#     f = []                  # Initialize Allele Frequency
#     for i in range(M):      # Iterate Though Alleles
#         f.append(af[an[i]]) # Get Frequency of Allele i
#     rho = compute_rho(f)    # Compute Individual Frequencies
    
#     # COMPUTE SIGMA FOR EACH ALLELE:
#     sg = []
#     for i in range(M):
#         ai = [j for j in range(MN) if a[j] == an[i]]
#         sg.append(compute_sigma([b[j] for j in ai],[p[j] for j in ai],es,pw))
#     I = compute_I(sg)
        
#     # COMPUTE SIGMA DOMINANT FOR EACH ALLELE:
#     sd = []
#     for i in range(M):
#         ai = [j for j in range(MN) if a[j] == an[i] & p[j] in idp]
#         sd.append(compute_sigma([b[j] for j in ai],[p[j] for j in ai],es,pw))
#     I = compute_I(sd)
    
    
        