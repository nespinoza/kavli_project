import disk_models
import numpy as np
import subprocess,os

Myr = 1e6

only_consider_these = 'e- H H2 He O C N Mg Si Fe \nS AL Ca Na Ni P K Ti CO\n'+\
                'OH SH N2 O2 SiO TiO SiS \nH2O C2 CH CN CS SiC NH SiH NO\n'+\
                'SN SiN SO S2 C2H HCN C2H2,acetylene \nCH4 ALH ALOH AL2O CaOH MgH MgOH PH3\n'+\
                'CO2 TiO2 Si2C SiO2 FeO NH2 NH3 CH2 CH3\nH2S KOH NaOH H+ H- Na+ K+ O2\n'+\
                'VO VO2 NaCL KCL FeH\nFe(L) MgSiO3(II) MgSiO3(III)\n'+\
                'AL2O3(L) KCL(L) SiC(L)\nNa2S(L) H2O(L)'

def checkCEA(compiler='gfortran'):
    """
    Check that all the fortran CEA codes have their 
    corresponding executables:
    """
    cwd = os.getcwd()
    os.chdir('CEA+FORTRAN')
    if not os.path.exists('FCEA2'):
        subprocess.Popen(compiler+' cea2.f -o FCEA2',\
                         shell = True).wait()
    if not os.path.exists('b1b2b3'):
        subprocess.Popen(compiler+' b1b2b3.f -o b1b2b3',\
                         shell = True).wait()
    if not os.path.exists('syntax'):
        subprocess.Popen(compiler+' syntax.f -o syntax',\
                         shell = True).wait()
    if not os.path.exists('thermo.lib'):
        p = subprocess.Popen(['./FCEA2'], stdout=subprocess.PIPE, \
                             stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_data = p.communicate(input='thermo')[0]
    if not os.path.exists('trans.lib'):
        p = subprocess.Popen(['./FCEA2'], stdout=subprocess.PIPE, \
                             stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_data = p.communicate(input='trans')[0]

    os.chdir(cwd)
    print '\t CEA up and running!'

def runCEA(filename):
    filename_no_extension = (filename.split('.inp')[0]).split('/')[1]
    subprocess.Popen('cp '+filename+' CEA+FORTRAN/.',shell = True).wait()
    cwd = os.getcwd()
    os.chdir('CEA+FORTRAN')
    p = subprocess.Popen(['./FCEA2'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_data = p.communicate(input=filename_no_extension)[0]
    subprocess.Popen('mv '+filename_no_extension+'.out ../CEAoutput/.',shell = True).wait() 
    subprocess.Popen('rm '+filename_no_extension+'.inp',shell = True).wait()
    os.chdir(cwd)

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
    
def calcCEA(T,P,name,N):
    """
    This function generates input files for CEA, runs them, and 
    spits out the equilibrium compositions:
    """

    f = open('CEAoutput/in_'+str(T)+'_'+str(P)+'.inp','w')
    f.write('# problem dataset: disk\n')
    f.write('problem tp ions\n')
    f.write('  p(bar) = '+str(P)+'\n')
    f.write('  t(k) = '+str(T)+'\n')
    f.write('reac\n')
    for i in range(len(name)):
        f.write('  na '+name[i]+' moles = '+str(N[i]).upper()+'\n')
    f.write('only '+only_consider_these)
    f.write('output calories siunits trace 1e-13\n')
    f.write('end')
    f.close()
    runCEA('CEAoutput/in_'+str(T)+'_'+str(P)+'.inp')

def readCEA(T,P):
    f = open('CEAoutput/in_'+str(T)+'_'+str(P)+'.out','r')
    species = []
    moles = []
    while True:
        line = f.readline()
        if 'MOLE FRACTIONS' in line:
            f.readline()
            break
        if line[:21] == ' CALCULATIONS STOPPED':
            f.close()
            return species,moles
    while True:
        line = f.readline()
        vec = line.split()
        if len(vec) == 0:
            break
        elif line[:21] == ' CALCULATIONS STOPPED':
            break
        else:
            species.append(vec[0])
            moles.append(np.double(vec[1]))
    f.close()
    return species,moles
checkCEA()
if not os.path.exists('CEAoutput'):
    os.mkdir('CEAoutput')
# Sample the disk coarsely in log:
r = 10**(np.linspace(-2,1.,100))
# Extract T-P profile and surface density:
T,P,Sigma = disk_models.ChambersModel(r,Myr)
# Extract abundances (note hydrogen is 1e12):
Z,name,N = read_abundances('mikeabundances.dat')
# Calculate equilibrium chemistry for each T-P point:
results = {}
for i in range(len(T)):
    calcCEA(T[i],P[i],name,N)
    c_species,c_moles = readCEA(T[i],P[i])
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
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.xlabel('Radius (AU)')
plt.ylabel('Mole abundance (normalized)')
plt.xscale('log')
plt.xlim([0.01,10.])
plt.ylim([1e-5,1.1])
for species in ['*N2','NH3','*CO','*CO2','H2O','*H','*H2']:
    plt.plot(r,results[species]/np.max(results[species]),\
                    linewidth=3,alpha=0.5,label=species)
plt.legend()
plt.show()         
#generate_inputCEA(T,P,name,logN,fname)

#import matplotlib.pyplot as plt
#plt.style.use('ggplot')
#plt.plot(Z,N,'o-')
#plt.xlabel('$Z$')
#plt.ylabel('$\log N(X)/N(H)$')
#plt.show()
#runCEA('example1.inp')
