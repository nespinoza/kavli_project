# -*- coding: utf-8 -*-
from pyx import *
import os
from scipy import interpolate
import disk_models
import utils
import numpy as np
import pickle 
Myr = 1e6

def emulate_colorbar_horizontal(g,teff,half_val,m,n,tick_range = None):
    # Define plot that will emulate the colorbar. First, define a fake y-axis:
    mymanualticks = [graph.axis.tick.tick(-2,label=""),graph.axis.tick.tick(-1,label="")]

    # Now a true y-axis with values I want to show:
    if tick_range is None:
        tick_range = [500,1200,1900,2600,3300]
        value_range = range(len(tick_range)-1)
    paint_tick_range = np.linspace(np.min(teff),np.max(teff),100)
    paint_value_range = np.double(np.arange(len(paint_tick_range)))
    paint_value_range = (paint_value_range/paint_value_range[-1])*3.0
    mymanualticks_x = []
    mymanualticks_x2 = []
    for i in range(len(value_range)):
            mymanualticks_x.append(graph.axis.tick.tick(value_range[i],label=str(tick_range[i])))
            mymanualticks_x2.append(graph.axis.tick.tick(value_range[i],label=""))

    cbar = c.insert(graph.graphxy(height=0.5,width=7, ypos = g.ypos + 7.5, \
                                                      xpos = g.xpos,\
           key=graph.key.key(pos=legend_pos,textattrs=[text.size(-3)]),\
           y = graph.axis.linear(min=-2,max=-1,title = None,manualticks=mymanualticks),\
           x2 = graph.axis.linear(min=np.min(value_range)-0.5,max=np.max(value_range)+0.5,\
                                  manualticks = mymanualticks_x,title = r'$T_\textnormal{disk}$',parter=None)))

    # "Paint" the plot:

    for i in range(len(paint_tick_range)):
       the_color = get_color(paint_tick_range,i,half_val,m,n)
       draw_rectangular_band(cbar, paint_value_range[i] - 0.5, -2, 1., 1.0, the_color)

    # Plot borders of cb
    cbar = c.insert(graph.graphxy(height=0.5,width=7, ypos = g.ypos + 7.5, \
                                                      xpos = g.xpos,\
           key=graph.key.key(pos=legend_pos,textattrs=[text.size(-3)]),\
           y = graph.axis.linear(min=-2,max=-1,title = None,manualticks=mymanualticks),\
           x2 = graph.axis.linear(min=np.min(value_range)-0.5,max=np.max(value_range)+0.5,\
                                  manualticks = mymanualticks_x2,title = None,parter=None)))

def emulate_colorbar_vertical(g,teff,half_val,m,n,tick_range = None):
    # Define plot that will emulate the colorbar. First, define a fake y-axis:
    mymanualticks = [graph.axis.tick.tick(-2,label=""),graph.axis.tick.tick(-1,label="")]

    # Now a true y-axis with values I want to show:
    if tick_range is None:
        tick_range = [500,1000,1500,2000,2500]
        value_range = range(len(tick_range)-1)
    paint_tick_range = np.linspace(np.min(teff),np.max(teff),100)
    paint_value_range = np.double(np.arange(len(paint_tick_range)))
    paint_value_range = (paint_value_range/paint_value_range[-1])*3.0-0.5
    mymanualticks_y = []
    mymanualticks_y2 = []
    for i in range(len(value_range)):
            mymanualticks_y.append(graph.axis.tick.tick(value_range[i],label=str(tick_range[i])))
            mymanualticks_y2.append(graph.axis.tick.tick(value_range[i],label=""))

    cbar = c.insert(graph.graphxy(height=5,width=0.5, ypos = g.ypos, \
                                                      xpos = g.xpos + 7.5,\
           key=graph.key.key(pos=legend_pos,textattrs=[text.size(-3)]),\
           x = graph.axis.linear(min=-2,max=-1,title = None,manualticks=mymanualticks),\
           y2 = graph.axis.linear(min=np.min(value_range)-0.5,max=np.max(value_range)+0.5,\
                                  manualticks = mymanualticks_y,title = r'$T_\textnormal{disk}$',parter=None)))

    # "Paint" the plot:

    for i in range(len(paint_tick_range)):
       the_color = get_color(paint_tick_range,i,half_val,m,n)
       draw_rectangular_band(cbar, -2, paint_value_range[i], 1., 1.0, the_color)

    # Plot borders of cb
    cbar = c.insert(graph.graphxy(height=5,width=0.5, ypos = g.ypos, \
                                                      xpos = g.xpos + 7.5,\
           key=graph.key.key(pos=legend_pos,textattrs=[text.size(-3)]),\
           x = graph.axis.linear(min=-2,max=-1,title = None,manualticks=mymanualticks),\
           y2 = graph.axis.linear(min=np.min(value_range)-0.5,max=np.max(value_range)+0.5,\
                                  manualticks = mymanualticks_y2,title = None,parter=None)))

def draw_rectangular_band(g, x_coord_init, y_coord_init, w, h, band_color):

    x0,y0 = g.pos(x_coord_init , y_coord_init)
    xf,yf = g.pos(x_coord_init+w, y_coord_init+h)

    x = x0
    y = y0
    w = np.abs(x0-xf)
    h = np.abs(y0-yf)

    g.fill(path.rect(x,y,w,h), [band_color])

def get_color(teff,i,half_val,m,n):
    # x goes from 1 to 0.5 to 0 from higher teff to lower teff:    
    x = (m*(teff[i])+n)    
    if teff[i] < half_val:
        c_color = x
        m_color = 0.5 + (x**2 - (3./2.)*x + 0.5)
        y_color = 1. - x
        k_color = 0.
    else:
        c_color = x
        m_color = 0.5 + (x**2 - (3./2.)*x + 0.5)
        y_color = 1. - x
        k_color = 0.

    print c_color,m_color,y_color,k_color
    return color.cmyk(c_color,m_color,y_color,k_color)

from pyx import text as pyx_text
from pyx import trafo
def draw_text(g, x_coord, y_coord, text_input, text_size = -2, color = None, rotation = 0.):
    """
    Function that draws text in a given plot
    INPUTS:
        g       (Object) A graph-type object to which you want to add the text
        x_coord     (Double) x-coordinate (in plot units) at which you want to place
                the text
        y_coord     (Double) y-coordinate (in plot units) at which you want to place
                the text
        text_input  (String) Text that you want to add.
        text_size   (int, optional) Text size of the text added to the plot. Default is -2.
        color       (instance) Color instance that defines the color that you want the 
                text to have. Default is black.
    """

    # First define the text attributes:
    textattrs = [pyx_text.size(text_size),pyx_text.halign.center, pyx_text.vshift.middlezero, trafo.rotate(rotation)]

    # Now convert plot positions to pyx's:
    x0,y0 = g.pos(x_coord, y_coord)

    # If no color is given, draw black text. If color is given, draw text with the input color:
    if color is None:
        g.text(x0,y0,text_input,textattrs)
    else:
        # First, check which was the input color palette:
        color_dict = color.color
        if len(color_dict.keys()) == 4:
            color_string = str(color_dict['c'])+','+str(color_dict['m'])+','+str(color_dict['y'])+','+str(color_dict['k'])
            color_palette = 'cmyk'
        else:
                        color_string = str(color_dict['r'])+','+str(color_dict['g'])+','+str(color_dict['b'])
                        color_palette = 'rgb'
        # Now draw the text:
        g.text(x0, y0, r"\textcolor["+color_palette+"]{"+color_string+"}{"+text_input+"}",textattrs)

# Extract data:
#radii,times,results = utils.get_data()
results = pickle.load(open('elemental_abundances_sim.pkl','r'))
radii = results['r']
times = results.keys()
times.remove('r')
times = np.double(np.array(times))
idx = np.argsort(times)
times = times[idx]
solar_abundances = pickle.load(open('elemental_abundances_ssn.pkl','r'))

# Plotting stuff:
unit.set(xscale = 0.8)
text.set(mode="latex")
text.preamble(r"\usepackage{color}")
text.preamble(r"\usepackage{wasysym}")
legend_text_size = -3
# More options on the legend:
min_x = np.min(radii)
max_x = np.max(radii)
legend_pos = 'ml'
xaxis = r'Radius (AU)'
yaxis = r'Fraction of element bound in gas/solids'

counter = 0
for i in range(len(times)):
    occupied = []
    c = canvas.canvas()
    if counter == 0:
        counter = 1
        min_y = -0.1
        max_y = 1.1

    print i,times[i]
    elemental_abundances = results[times[i]]

    g = c.insert(graph.graphxy(height=5,width=7, key=graph.key.key(pos=legend_pos,\
                                             textattrs=[text.size(legend_text_size)]),\
                   x = graph.axis.logarithmic(min=min_x,max=max_x,title = xaxis, texter=graph.axis.texter.decimal()),\
                   y = graph.axis.linear(min=min_y,max=max_y,title = yaxis)))


    draw_text(g, radii[int(len(radii)/2.)], max_y + max_y*0.1, '{0: 2.2f} Myr disk'.format(times[i]), text_size = 2)
    for element in solar_abundances.keys():
        gas = elemental_abundances[element+'(g)']/solar_abundances[element]
        solid = elemental_abundances[element+'(s)']/solar_abundances[element]
        if element in ['Na','K','Ti','V','O','Fe']:
                if element == 'Fe':
                   clr = color.cmyk.RoyalBlue
                elif element == 'O':
                   clr = color.cmyk.Thistle
                elif element == 'V':
                   clr = color.cmyk.Black
                elif element == 'Ti':
                   clr = color.cmyk.Gray
                elif element == 'K':
                   clr = color.cmyk.OliveGreen
                elif element == 'Na':
                   clr = color.cmyk.Orange
                color_dict = clr.color
                color_string = str(color_dict['c'])+','+str(color_dict['m'])+','+str(color_dict['y'])+','+str(color_dict['k'])
                color_palette = 'cmyk'
                g.plot(graph.data.values(x=radii, y=solid, title = r"\textcolor["+color_palette+"]{"+color_string+"}{"+element+"}"),\
                                         styles = [graph.style.line([clr,\
                                         style.linestyle.solid,\
                                         style.linewidth.thick])])

                g.plot(graph.data.values(x=radii, y=gas, title = None),\
                                         styles = [graph.style.line([clr,\
                                         style.linestyle.dashed,\
                                         style.linewidth.thick])])
    if i<10:
        c.writeGSfile("colormap_00"+str(i)+".png",resolution=1000)
    elif i<100:
        c.writeGSfile("colormap_0"+str(i)+".png",resolution=1000)
    else:
        c.writeGSfile("colormap_"+str(i)+".png",resolution=1000)

