# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pyCEA
import disk_models
import numpy as np
import subprocess,os

Myr = 1e6

only_consider_these =  '\n  e- AL AL2O ALH ALO ALOH ALS\n'+\
                       '  C CH4 CN CO CO2 CP CS\n'+\
                       '  Ca CaH CaO CaOH CaS\n'+\
                       '  Cr\n'+\
                       '  Fe FeO\n'+\
                       '  H H2 H2O H2S HCN HCO HPO\n'+\
                       '  He\n'+\
                       '  Mg MgH MgN MgO MgOH MgS\n'+\
                       '  N NH3 NO\n'+\
                       '  Na Na2 NaH NaO NaOH\n'+\
                       '  Ni NiO NiS\n'+\
                       '  O O2\n'+\
                       '  P PH PN PO PS\n'+\
                       '  S S2 SN SO SO2\n'+\
                       '  Si SiC SiH SiN SiO SiS\n'+\
                       '  Ti TiO TiO2\n'+\
                       '  AL2O3(a)  ALN(cr)\n'+\
                       '  C(gr)\n'+\
                       '  CaS(cr)\n'+\
                       '  Cr(cr)\n'+\
                       '  Fe(a) Fe(c) Fe(d)\n'+\
                       '  H2O(cr)\n'+\
                       '  MgAL2O4(cr) MgS(cr) Mg2SiO4(cr) MgSiO3(I) MgSiO3(II) MgSiO3(III)\n'+\
                       '  Ni(cr) P(cr) Si(cr) SiC(b)\n'+\
                       '  TiC(cr) TiN(cr) Ti2O3(I)\n'

# Times to sample:
all_times = np.linspace(0.01,3.,100)

# Sample the disk coarsely in log:
r = 10**(np.linspace(-2,0.1,50))

# Extract abundances (note hydrogen is 1e12):
Z,name,N = pyCEA.read_abundances('ssabundances_4Gyr.dat')

# Chech which elements we are actually going to use:
idx = pyCEA.check_elements(name,only_consider_these)
Z = Z[idx]
name = name[idx]
N = N[idx]

for m in range(len(all_times)):
  t = all_times[m]
  # Extract T-P profile and surface density:
  T,P,Sigma = disk_models.ChambersModel(r,t*Myr)
  # Calculate equilibrium chemistry for each T-P point:
  results = {}
  for i in range(len(T)):
      pyCEA.calcCEA(T[i],P[i],name,N,only_consider_these,prefix='benchmark')
      c_species,c_moles = pyCEA.readCEA(T[i],P[i],prefix='benchmark')
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

  plt.style.use('ggplot')
  plt.title('{0:2.2} Myr Disk'.format(t))
  plt.xlabel('Radius (AU)')
  plt.ylabel('Relative molar abundance')
  plt.xscale('log')
  plt.xlim([0.01,0.4])
  plt.ylim([0.0,1.8e-7])
  print results.keys()
  for species in ['*Ti','*TiO','TiO2','TiN(cr)']:
      plt.plot(r,results[species],\
                      linewidth=3,alpha=0.5,label=species)
  plt.legend(loc=4)
  if m<10:
      plt.savefig('disk_Ti_00'+str(m)+'.png')
  elif m<100:
      plt.savefig('disk_Ti_0'+str(m)+'.png')
  else:
      plt.savefig('disk_Ti_'+str(m)+'.png')
  plt.clf()
  #plt.show()         
