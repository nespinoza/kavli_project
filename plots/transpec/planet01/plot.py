from scipy.ndimage.filters import gaussian_filter
import numpy as np
from pyx import *

# Plotting stuff:
unit.set(xscale = 0.8)
text.set(mode="latex")
text.preamble(r"\usepackage{color}")
text.preamble(r"\usepackage{wasysym}")
legend_text_size = -3

wssn,rssn = np.loadtxt('ssn_acc_0.1_1X.dat',unpack=True,usecols=(0,1))
w,r = np.loadtxt('acc_0.1_1X.dat',unpack=True,usecols=(0,1))

# More options on the legend:
min_y = None#np.min(r)
max_y = 0.1295#np.max(r)
min_x = 0.5
max_x = 2.0
legend_pos = 'tr'
Rstar = 6.957e5 # km, solar radii
xaxis = r'Wavelength (microns)'
yaxis = r'$R_p/R_*$'
idx = np.where((w>min_x)&(w<max_x))[0]
c = canvas.canvas()
g = c.insert(graph.graphxy(height=5,width=10, key=graph.key.key(pos=legend_pos,\
                                             textattrs=[text.size(legend_text_size)]),\
                   x = graph.axis.linear(min=min_x,max=max_x,title = xaxis, texter=graph.axis.texter.decimal()),\
                   y = graph.axis.linear(min=min_y,max=max_y,title = yaxis)))

g.plot(graph.data.values(x=wssn[idx], y=gaussian_filter((rssn/Rstar)[idx],3), title = 'SSN abundances'),\
                             styles = [graph.style.line([color.cmyk.Gray,\
                                       style.linestyle.solid,\
                                       style.linewidth.thin])])

g.plot(graph.data.values(x=w[idx], y=gaussian_filter((r/Rstar)[idx],3), title = 'Accreted abundance'),\
                             styles = [graph.style.line([color.cmyk.RoyalBlue,\
                                       style.linestyle.solid,\
                                       style.linewidth.thin])])

print r/Rstar
c.writeGSfile("transpec.png",resolution=1000)
