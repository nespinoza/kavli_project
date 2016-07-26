# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
import utils
import pickle

radius,times,results = utils.get_data()
Z,name,N = utils.read_abundances('ssabundances_4Gyr.dat')
abundances = {}
for i in range(len(name)):
    abundances[name[i]] = N[i]
Habundance = abundances['H']
for i in range(len(name)):
    abundances[name[i]] = abundances[name[i]]/Habundance
element1 = 'Na'
element2 = 'K'
print 'Initial ratio:',abundances[element1]/abundances[element2]
psolar = abundances[element1]/abundances[element2]
idx_time = 0.5
for i in range(len(times)):
  i = int(len(times)/2.)
  for element1 in abundances.keys():
    t = times[i]
    resulting_abundances = results[t]
    solids_element1 = 0.0
    gas_element1 = 0.0
    solids_element2 = 0.0
    gas_element2 = 0.0
    for species in resulting_abundances.keys():
        if element1 in species:
            # Check how many atoms of the element we have in the species:    
            fac = utils.howmany(element1,species)
            # Check if solid or gas:
            if '(a)' in species or '(cr)' in species or '(b)' in species\
               or '(c)' in species or '(I)' in species or "(I')" in species or\
               'II' in species:
                solids_element1 = solids_element1 + fac*resulting_abundances[species]
            else:
                gas_element1 = gas_element1 + fac*resulting_abundances[species]

    total_element1 = gas_element1 + solids_element1
    # Transform to percentage of gas and solids (in number):
    if element1 in ['Na','K','Ti','V','O','N','Fe']:
        lin = plt.plot(radius,gas_element1/total_element1,'-',linewidth = 2,label = element1)
        plt.plot(radius,solids_element1/total_element1,'--',linewidth = 2,color = lin[0].get_color())
    #gas_element1 = (gas_element1/total_element1)*abundances[element1]
    #solids_element1 = (solids_element1/total_element1)*abundances[element1] 
  break
plt.title('5 Myr disk')
plt.xlabel('Radius (AU)')
plt.ylabel(r'Elemental abundance ($\times$ solar)')
plt.ylim([-0.1,1.1])
plt.xscale('log')
plt.legend()
plt.show()
