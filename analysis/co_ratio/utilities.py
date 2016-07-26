import periodictable as pt
import numpy as np


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

def realZ():
    """
    Returns the real Z given the Lodders compilation
    """

    AN,name,N = read_abundances('ssabundances_4Gyr.dat')
    num = 0.0
    den = 0.0
    for i in range(len(name)):
        element = name[i]
        # Calculate metals:
        if element != 'H' and element != 'He':
            if element == 'CL':
                exec 'num = num + (N[i])*pt.Cl.mass'
            elif element == 'AL':
                exec 'num = num + (N[i])*pt.Al.mass'
            else:
                exec 'num = num + (N[i])*pt.'+element+'.mass'
        if element == 'CL':
            exec 'den = den + (N[i])*pt.Cl.mass'
        elif element == 'AL':
            exec 'den = den + (N[i])*pt.Al.mass'
        else:
            exec 'den = den + (N[i])*pt.'+element+'.mass'
    return num/den
def getZ(ratio, solids=False, gas = False, no = 1.):
    """
    Given a C/O ratio, this returns the metallicity 
    of the compund.
    """
    AN,name,N = read_abundances('ssabundances_4Gyr.dat')
    for i in range(len(name)):
        if name[i] == 'O':
            NO = N[i]
            break
    C1 = 0.0
    C2 = 0.0
    for i in range(len(name)):
        element = name[i]
        # Calculate C1 and C2:
        if element != 'H' and element != 'He' and \
           element != 'C' and element != 'O':
            if element == 'CL':
                exec 'C1 = C1 + (N[i]/NO)*pt.Cl.mass'
            elif element == 'AL':
                exec 'C1 = C1 + (N[i]/NO)*pt.Al.mass'
            else:
                exec 'C1 = C1 + (N[i]/NO)*pt.'+element+'.mass'
        if element != 'C' and element != 'O':
            if element == 'CL':
                exec 'C2 = C2 + (N[i]/NO)*pt.Cl.mass'
            elif element == 'AL':
                exec 'C2 = C2 + (N[i]/NO)*pt.Al.mass'
            else:
                if solids:
                    if element != 'H' and element != 'He':
                        exec 'C2 = C2 + (N[i]/NO)*pt.'+element+'.mass'
                else:
                    exec 'C2 = C2 + (N[i]/NO)*pt.'+element+'.mass'
                #exec 'C2 = C2 + (N[i]/NO)*pt.'+element+'.mass'
    if gas:
        C1 = 0.0
        C2 = pt.H.mass/(no*1e-4)
       
    num = ratio*pt.C.mass + pt.O.mass + C1
    den = ratio*pt.C.mass + pt.O.mass + C2
    #print num
    #print den
    return num/den
