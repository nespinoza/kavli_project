# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pyCEA
import disk_models
import numpy as np
import subprocess,os

Myr = 1e6

only_consider_these =  '\n  AL ALN ALC ALH ALH2 ALO ALOH ALS\n'+\
                       '  C CH4 CN CO CO2 CS\n'+\
                       '  Ca CaH CaO CaOH CaS\n'+\
                       '  K KCN KH KOH\n'+\
                       '  Fe FeH FeO FeS Fe(OH)2\n'+\
                       '  H H2 H2O H2S HCN HCO HS\n'+\
                       '  He\n'+\
                       '  Mg MgH MgN MgO MgOH MgS\n'+\
                       '  N N2 NH NH3 NO\n'+\
                       '  Na Na2 NaH NaO NaOH\n'+\
                       '  O O2 OH \n'+\
                       '  S S2 SN SO SO2\n'+\
                       '  Si SiH SiH4 SiN SiO SiS\n'+\
                       '  Ti TiO TiO2 TiS\n'+\
                       '  V VN VO VO2 \n'+\
                       '  NaALSi3O8(cr)\n'+\
                       '  MgAL2O4(cr)\n'+\
                       '  AL2O3(a)\n'+\
                       '  Ca2AL2SiO7(cr) Mg2SiO4(cr)\n'+\
                       '  FeS(a) FeS(b) FeS(c)\n'+\
                       '  KALSi3O8(cr)\n'+\
                       '  CaAL12O19(cr) CaAL2Si2O8(cr)\n'+\
                       '  CaTiO3(a)\n'+\
                       '  V2O3(cr)\n'+\
                       '  Fe(a) Fe(c) Fe3O4(cr)\n'+\
                       "  TiN(cr) Ti2O3(I) Ti2O3(I')\n"+\
                       '  MgS(cr) MgSiO3(II) AL(cr) Ca2AL2SiO7(cr)\n'
#                       "  TiN(cr) Ti2O3(I) Ti2O3(I')\n"+\
#                       '  MgS(cr) MgSiO3(II) AL(cr) Ca2AL2SiO7(cr)\n'
#                       '  ALN(cr)\n'+\
#                       '  MgS(cr) MgSiO3(II)\n'
#                       '  CaAL12O19(cr) CaAL2Si2O8(cr)\n'+\                       
#                       '  CaAL2Si2O8(cr) ALN(cr)\n'+\
#                       '  MgS(cr) MgSiO3(II)\n'
#                       '  Si(cr) SiC(b) \n'+\
#                       '  Na8AL6Si6O24CL2(cr)\n'
#                       ''
#                      '  MgS(cr) MgSiO3(I) MgSiO3(II) MgSiO3(III)\n'+\
#                       '  CaS(cr)\n'+\
#                       '  Ca3P2O8(cr) Ca2AL2SiO7(cr)\n'
#                       '  CaMgSi2O6(cr)\n'+\
#                       '  MgAL2O4(cr) MgS(cr) Mg2SiO4(cr) MgSiO3(I) MgSiO3(II) MgSiO3(III)\n'+\
#                       '  NaALSi3O8(cr) \n'+\
#                       '  Fe3C(cr) Fe3P(cr) Fe2SiO4(a) Fe2SiO4(b)\n'+\
#                       '  FeCr2O4(a) FeCr2O4(b)\n'+\
#                       '  FeS(a) FeS(b) FeS(c) FeSiO3(a) FeSiO3(b)\n'
#                       '  H2O(cr)\n'+\
#                       '  KALSi3O8(cr) KALSi3O8(L)\n'+\

def erase_from_list(sp,long_string):
    vec_elements = long_string.split()
    if '*' in sp:
        sp = sp[1:]
    vec_elements.remove(sp)
    out_string = '\n  '
    idx_space = 0 
    for i in range(len(vec_elements)):
        out_string = out_string+vec_elements[i]
        idx_space = idx_space + 1
        if idx_space == 5:
            out_string = out_string+'\n  '
            idx_space = 0
        else:
            if i == len(vec_elements)-1:
                out_string = out_string+'\n'
            else:
                out_string = out_string+' '
    return out_string

# Times to sample:
all_times = np.linspace(0.01,10.,100)
#all_times = np.linspace(0.01,10.,100)

# Sample the disk coarsely in log:
#r = 10**(np.linspace(-2,-0.8,500))
r = 10**(np.linspace(-2,-0.8,100))

# Extract abundances (note hydrogen is 1e12):
Zo,nameo,No = pyCEA.read_abundances('ssabundances_4Gyr.dat')

# Chech which elements we are actually going to use:
idx = pyCEA.check_elements(nameo,only_consider_these)
Z = Zo[idx]
name = nameo[idx]
N = No[idx]
all_results = {}
all_results['r'] = r
for m in range(len(all_times)):
  print 'Time stamp ',m,' of ',len(all_times),'...'
  non_zero_species = []
  t = all_times[m]
  # Extract T-P profile and surface density:
  T,P,Sigma = disk_models.ChambersModel(r,t*Myr)
  P = np.ones(len(P))*1e-4
  # Calculate equilibrium chemistry for each T-P point:
  results = {}
  only_consider_these_now = only_consider_these
  for i in range(len(T)):
      pyCEA.calcCEA(T[i],P[i],name,N,only_consider_these_now,prefix='final')
      c_species,c_moles = pyCEA.readCEA(T[i],P[i],prefix='final')
      if i == 0:
          for j in range(len(c_species)):
              results[c_species[j]] = np.zeros(len(r))
              results[c_species[j]][i] = c_moles[j]
      else:
          for j in range(len(c_species)):
              if c_species[j] in results.keys():
                  results[c_species[j]][i] = c_moles[j]
              else:
                  results[c_species[j]] = np.zeros(len(r))
                  results[c_species[j]][i] = c_moles[j]
  all_results[t] = results

FILE = open('results.pkl','w')
import pickle
pickle.dump(all_results, FILE)
FILE.close()
#  plt.legend(loc=4)
#  plt.show()
