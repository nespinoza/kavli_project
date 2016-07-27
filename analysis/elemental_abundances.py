# -*- coding: utf-8 -*-
import periodictable as pt
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
import utils
import pickle

radius,times,results = utils.get_data()
Z,name,N = utils.read_abundances('ssabundances_4Gyr.dat')
##### Sanity check #####
#for i in range(len(times)):
#  resulting_abundances = results[times[i]]
#  in_element = {}
#  out_element = {}
#  counter = 0
#  tot = 0.0
#  tot_abundances = 0.0
#  all_s = 0.0
#  idx_radius = 0
#  for species in resulting_abundances.keys():
#    all_s = all_s + resulting_abundances[species][idx_radius]
#    for element in ['H','He','C','N','O','Na','Mg','AL','Si','S','K','Ca','Ti','V','Fe']:
#        if counter == 0:
#            idx = np.where(name == element)[0]
#            tot_abundances = tot_abundances + N[idx]
#            in_element[element] = N[idx]
#            out_element[element] = 0.0
#        fac = utils.howmany(element,species)
#        out_element[element] = out_element[element] + fac*resulting_abundances[species][idx_radius]
#        tot = tot+fac*resulting_abundances[species][idx_radius]
#        if element == 'Fe':
#            print species,resulting_abundances[species][idx_radius],'(',fac,')'
#    counter = 1
#  print 'all_s:'
#  print all_s
#  print 'tot:'
#  print tot
#  totcheck = 0.0
#  for element in ['H','He','C','N','O','Na','Mg','AL','Si','S','K','Ca','Ti','V','Fe']:
#    totcheck = totcheck + out_element[element]
#    print element+' in:  ',in_element[element]/tot_abundances
#    print element+' out: ',out_element[element]/tot
#  break
#sys.exit()
########################
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
