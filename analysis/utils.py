# -*- coding: utf-8 -*-
from scipy import interpolate
import pickle
import numpy as np

def get_accretion_data(filename):
    fin = open(filename,'r')
    time = np.array([])
    a = np.array([])
    Mg = np.array([])
    TM = np.array([])
    a_width = np.array([])
    while True:
        line = fin.readline()
        if line != '':
            if line[0] != '#':
                vec = line.split()
                time = np.append(time,np.double(vec[0].replace('D','e')))
                a = np.append(a,np.double(vec[1].replace('D','e')))
                Mg = np.append(Mg,np.double(vec[2].replace('D','e')))
                TM = np.append(TM,np.double(vec[3].replace('D','e')))
                a_width = np.append(a_width,np.double(vec[4].replace('D','e')))
        else:
            break
    fin.close()
    return time,a,Mg,TM,a_width

def read_abundances(filename):
    """
    This code reads abundances from a file containing Z, name, A(X) and error on this 
    number, where A(X) = log N(X)/N(H) + 12, and N(X) is the number of atoms of element 
    X. The code returns:

       Z            Atomic number of the element X
       Name         Name of the element X
       10**A        10 to the (log N(X)/N(H) - 12), where N(X) is the number of atoms of element X,
                    and N(H) is the number of hydrogen atoms (thus, this takes 1e12 for hydrogen 
                    by definition).
    """
    Z = np.array([])
    name = np.array([])
    logN = np.array([])

    f = open(filename,'r')
    while True:
        line = f.readline()
        if line != '':
            if line[0] != '#':
                data = line.split()
                Z = np.append(Z,np.double(data[0]))
                name = np.append(name,data[1])
                logN = np.append(logN,np.double(data[2]))
        else:
            break
    f.close()
    return Z,name,10**(logN)

def checkelement(element, species):
    if element == 'I':
        return False
    idx = species.find(element)
    if idx != -1:
        l = len(element)
        ls = len(species)
        if l == ls:
            return True
        else:
            s = species[idx:]
            if len(s) == l:
                return True
            if element == 'CL' or element == 'AL':
                if s[2].isdigit():
                    return True
                elif s[2].lower() == s[2]:
                    if len(s) == 3:
                        return False
                    else:
                        return checkelement(element,s[3:])
                else:
                    return True
            else:
                if s[l].upper() == s[l]:
                    if element == 'C' and s == 'CL':
                        return False
                    else: 
                        return True
                else:
                    return checkelement(element,s[l:])

    else:
        return False

def howmany(element,species):
    """
    Given an input element and a species, it returns the 
    number of elements in the species. For example, for 
    the species 
 
                            CH4

    It returns 1 if element is 'C', and 4 if element is 'H'.
    """
    if species == 'Fe(OH)2':
        species = 'FeO2H2'
        
    l = len(element)
    ls = len(species)
    if species.find(element) == -1:
        return 0.
    elif ls == l:
        return 1.
    else:
        for i in range(ls):
            if i+l == ls + 1:
               return 0.0
            if species[i:i+l] == element:
                if i+l == ls:
                    return 1.
                elif species[i+l:i+l+1].isdigit():
                    kicked_out = False
                    for k in range(1,ls-(i+l)+1):
                        if not species[i+l:i+l+k].isdigit():
                            kicked_out = True
                            break
                    if kicked_out:
                        return np.double(species[i+l:i+l+k-1])
                    else:
                        return np.double(species[i+l:i+l+k])
                else:
                    if species[i+l:i+l+1] == species[i+l:i+l+1].upper():
                        return 1.
    return 0.0

def read_abundances(filename):
    """
    This code reads abundances from a file containing Z, name, A(X) and error on this 
    number, where A(X) = log N(X)/N(H) + 12, and N(X) is the number of atoms of element 
    X. The code returns:

       Z            Atomic number of the element X
       Name         Name of the element X
       10**A        10 to the (log N(X)/N(H) - 12), where N(X) is the number of atoms of element X,
                    and N(H) is the number of hydrogen atoms (thus, this takes 1e12 for hydrogen 
                    by definition).
    """
    Z = np.array([])
    name = np.array([])
    logN = np.array([])

    f = open(filename,'r')
    while True:
        line = f.readline()
        if line != '':
            if line[0] != '#':
                data = line.split()
                Z = np.append(Z,np.double(data[0]))
                name = np.append(name,data[1])
                logN = np.append(logN,np.double(data[2]))
        else:
            break
    f.close()
    return Z,name,10**(logN)

def get_data():
    results = pickle.load(open('results.pkl','r'))
    # Extract data in user-friendly format:
    radius = results['r']
    times = results.keys()
    times.remove('r')
    times = np.array(times)
    idx = np.argsort(times)
    times = times[idx]
    # Post-processing in case of nans and zeros on all species
    # (i.e., failed calculations from CEA):
    for i in range(len(times)):
        abundances = results[times[i]]
        # First check for nans and sum abundances at the 
        # same time (if zero in all abundances, it means 
        # calculations failed for CEA):
        tot = 0.0
        for species in abundances.keys():
            idx = np.isnan(abundances[species])
            if len(idx) != len(abundances[species]):
                f = interpolate.interp1d(radius[~idx], abundances[species][~idx])
                abundances[species][idx] = f(radius[idx])
            tot = tot + abundances[species]
        idx = np.where(tot == 0.0)[0]
        if len(idx)>0:
            idx_non_zero = np.where(tot != 0.0)[0]
            for species in abundances.keys():
                f = interpolate.interp1d(radius[idx_non_zero], abundances[species][idx_non_zero])
                abundances[species][idx] = f(radius[idx])
    return radius,times,results
