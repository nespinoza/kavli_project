import matplotlib.pyplot as plt
import numpy as np

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

# Extract data from table:
data = read_latex_table()
# Get name, planet mass, radius, stellar metallicity and metal mass:
names = data[:,1]
mass = data[:,2]
radius = data[:,3]
metallicity = data[:,6]
metals = data[:,7]

print mass[0]
print read_latex_errors(mass[0])
