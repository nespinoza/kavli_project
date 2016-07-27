# -*- coding: utf-8 -*-
import periodictable as pt
import matplotlib.pyplot as plt
import utilities
import numpy as np

Mearth = 5.972e24 # kg
Mjup = 317.8*Mearth
amu = 1.660539040e-27 # kg

def read_latex_table(fname='data.tex'):
    f = open(fname,'r')
    row = ''
    first_time = True
    while True:
        line = f.readline()
        if line != '':
            if r'\\' in line:
                vec = line.split(r'\\')
                row = row + vec[0]
                # Got to final element of the row, separate:
                elements = row.split('&')
                if first_time:
                    output = np.array(elements)
                    first_time = False
                else:
                    output = np.vstack((output,np.array(elements)))
                row = ''
            else:
                row = row + line
        else:
            break
    return output

def read_latex_errors(s):
    s = s.split('$')[1]
    if r'\\pm' in s:
        val,err = s.split(r'\\pm')
        return np.double(val),np.double(err),np.double(err)
    elif r'\pm' in s:
        val,err = s.split(r'\pm')
        return np.double(val),np.double(err),np.double(err)
    else:
        if s.find('_')>s.find('^'):
            val,errors = s.split('^')
            err_up,err_down = errors.split('_')
        else:
            val,errors = s.split('_')
            err_down,err_up = errors.split('^')
        if '{' in err_up:
            err_up = err_up[1:-1]
            err_down = err_down[1:-1]
        return np.double(val),np.double(err_down),np.double(err_up)

# Estimate the constants in the expression to obtain total number of oxygen atoms:
Z,name,N = utilities.read_abundances('ssabundances_4Gyr.dat')
# Extract data from table:
data = read_latex_table()
# Get name, planet mass, radius, stellar metallicity and metal mass:
names = data[:,1]
mass = data[:,2]
radius = data[:,3]
metallicity = data[:,6]
metals = data[:,7]

val,err_down,err_up = read_latex_errors(mass[0])
r = np.linspace(1.,100.,10000)
T = 200.*((r)**(-0.62))
O_gas = np.zeros(len(r))
CO_ratio_gas = np.zeros(len(r))
CO_ratio_solids = np.zeros(len(r))
bulk_co = np.zeros(len(r))

print names[0]
Mcore = 10.
Msolids,errdown,errup = read_latex_errors(metals[0]) 
Msolids = Mearth*(Msolids-Mcore)
totalmass,errdown,errup = read_latex_errors(mass[0]) 
totalmass = Mjup*totalmass
Mgas = totalmass - Msolids
mo = pt.O.mass*amu
mc = pt.C.mass*amu
mh = pt.H.mass*amu
mhe = pt.He.mass*amu
for i in range(len(T)):
    temp = T[i]
    Csolids = 0.0
    Osolids = 0.0
    Cgas = 0.0
    Ogas = 0.0
    if temp<=20:
        Csolids = Csolids + 1.5
        Osolids = Osolids + 1.5
    else:
        Cgas = Cgas + 1.5
        Ogas = Ogas + 1.5
    if temp<=47:
        Csolids = Csolids + 0.3
        Osolids = Osolids + 0.6
    else:
        Cgas = Cgas + 0.3
        Ogas = Ogas + 0.6
    if temp<=135:
        Osolids = Osolids + 0.9
    else:
        Ogas = Ogas + 0.9
    if temp<=500:
        Csolids = Csolids + 0.6
    else:
        Cgas = Cgas + 0.6
    if temp<=1500:
        Osolids = Osolids + 1.4
    else:
        Ogas = Ogas + 1.4
    if Ogas != 0:
        O_gas[i] = Ogas
        CO_ratio_gas[i] = Cgas/Ogas

        CO_ratio_solids[i] = Csolids/Osolids
        Nogas = (Mgas/(mo + (CO_ratio_gas[i]*mc) + (1e4/Ogas)*(mh + 0.096827*mhe)))
        Nosolids = (Msolids)/(mo + mc*CO_ratio_solids[i] + 2.4e-26)
        No = Nogas + Nosolids
        bulk_co[i] = (1./mc)*(((Mgas+Msolids)/No) - mo - 3.87e-24)
print Msolids/Mearth
print totalmass/Mjup
from pylab import *
plot(r,bulk_co)
show()

