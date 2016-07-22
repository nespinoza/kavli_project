# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
import utils
import pickle

def get_data():
    results = pickle.load(open('results.pkl','r'))
    # Extract data in user-friendly format:
    radius = results['r']
    times = results.keys()
    times.remove('r')
    times = np.array(times)
    idx = np.argsort(times)
    times = times[idx]
    return radius,times,results

radius,times,results = get_data()
Z,name,N = utils.read_abundances('ssabundances_4Gyr.dat')
abundances = {}
for i in range(len(name)):
    abundances[name[i]] = N[i]
Habundance = abundances['H']
for i in range(len(name)):
    abundances[name[i]] = abundances[name[i]]/Habundance
element1 = 'Fe'
element2 = 'O'
print 'Initial ratio:',abundances[element1]/abundances[element2]
psolar = abundances[element1]/abundances[element2]
idx_time = 0.5
for i in range(len(times)):
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
        if element2 in species:
            # Check how many atoms of the element we have in the species:    
            fac = utils.howmany(element2,species)
            # Check if solid or gas:
            if '(a)' in species or '(cr)' in species or '(b)' in species\
               or '(c)' in species or '(I)' in species or "(I')" in species or\
               'II' in species:
                solids_element2 = solids_element2 + fac*resulting_abundances[species]
            else:
                gas_element2 = gas_element2 + fac*resulting_abundances[species]
    total_element1 = gas_element1 + solids_element1
    total_element2 = gas_element2 + solids_element2
    # Transform to percentage of gas and solids (in number):
    gas_element1 = (gas_element1/total_element1)*abundances[element1]
    solids_element1 = (solids_element1/total_element1)*abundances[element1] 
    gas_element2 = (gas_element2/total_element2)*abundances[element2]
    solids_element2 = (solids_element2/total_element2)*abundances[element2]
    # Check for zeroes in gas and solid elements:
    idx_gas = np.where(gas_element2 != 0.0)[0]
    gas_ratio = gas_element1[idx_gas]/gas_element2[idx_gas]
    idx_solid = np.where(solids_element2 != 0.0)[0]
    solid_ratio = solids_element1[idx_solid]/solids_element2[idx_solid]
    if t > idx_time:
        plt.plot(radius[idx_gas],gas_ratio/psolar,'--',label = str(t)+' Myr disk (gas)')
        plt.plot(radius[idx_solid],solid_ratio/psolar,'-',label = str(t)+' Myr disk (solid)')
        idx_time = idx_time + 2.0

plt.xlabel('Radius (AU)')
plt.ylabel(element1+'/'+element2+r' ratio ($\times$ solar)')
plt.xscale('log')
plt.legend()
plt.show()
