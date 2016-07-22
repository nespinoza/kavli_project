# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
import utils
import pickle

if not os.path.exists('elemental_abundances_sim.pkl'):
    # First, generate or load the solar-system nebula 
    # abundances (normalized to the number of hydrogen atoms):
    if not os.path.exists('elemental_abundances_ssn.pkl'):
        abundances = {}
        Z,name,N = utils.read_abundances('ssabundances_4Gyr.dat')
        # Normalize number abundances with respect to 
        # hydrogen:
        for i in range(len(name)):
            abundances[name[i]] = N[i]
        Habundance = abundances['H']
        for i in range(len(name)):
            abundances[name[i]] = abundances[name[i]]/Habundance
        # Save:
        FILE = open('elemental_abundances_ssn.pkl','w')
        pickle.dump(abundances, FILE)
        FILE.close()

    else:
        FILE = open('elemental_abundances_ssn.pkl','r')
        abundances = pickle.load(FILE)
        FILE.close()
        name = abundances.keys()

    # Get data from the simulations:
    radius,times,results = utils.get_data()
    # Define dictionary to save the elemental abundances:
    elemental_abundances = {}

    # Generate elemental abundances (by number,
    # in units of hydrogen abundance) dictionary; 
    # this will save, at a given time, the abundance 
    # in gas and solids for each element at each radius:
    for i in range(len(times)):
        elemental_abundances[times[i]] = {}
        for element in name:
            elemental_abundances[times[i]][element+'(g)'] = np.zeros(len(radius))
            elemental_abundances[times[i]][element+'(s)'] = np.zeros(len(radius))
    # Now, at each time, extract elemental abundances in gas and 
    # solids:    
    for i in range(len(times)):
        t = times[i]
        for element in name:
            resulting_abundances = results[t]
            solids = 0.0
            gas = 0.0
            for species in resulting_abundances.keys():
                if utils.checkelement(element, species):
                    # Check how many atoms of the element we have in the species:    
                    fac = utils.howmany(element,species)
                    # Check if solid or gas:
                    if '(a)' in species or '(cr)' in species or '(b)' in species\
                       or '(c)' in species or '(I)' in species or "(I')" in species or\
                       'II' in species:
                        solids = solids + fac*resulting_abundances[species]
                    else:
                        gas = gas + fac*resulting_abundances[species]
            total_element = gas + solids
            # Transform to percentage of gas and solids (in number, normalized to H 
            # number), and save. If total_element is zero (i.e., element not considered in 
            # our model), then give it all to the gas:
            if type(total_element) is not float:
                idx_zeros = np.where(total_element==0)[0]
                if len(idx_zeros)>0:
                    print element
                elemental_abundances[t][element+'(g)'] = (gas/total_element)*abundances[element]
                elemental_abundances[t][element+'(s)'] = (solids/total_element)*abundances[element]
            else:
                elemental_abundances[t][element+'(g)'] = abundances[element]

    FILE = open('elemental_abundances_sim.pkl','w')
    pickle.dump(elemental_abundances, FILE)
    FILE.close()
else:
    FILE = open('elemental_abundances_sim.pkl','r')
    elemental_abundances = pickle.load(FILE)
    FILE.close()

# Now load planet: 


