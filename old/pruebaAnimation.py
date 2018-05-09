# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 01:52:15 2018

@author: Stradi
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import numpy as np # just for dummy data generation

# generate dummy data
ndat = 20
x = np.linspace(0,1,ndat)
phi = np.linspace(0,2*np.pi,100,endpoint=False)
dat = np.transpose([x[:,None]*np.cos(phi),x[:,None]*np.sin(phi)],(1,2,0))

# create figure and axes
fig = plt.figure()
ax_pl = plt.subplot2grid((5,5),(0,0),colspan=5,rowspan=3)  # axes_plot
ax_bl = plt.subplot2grid((5,5),(4,0),colspan=2,rowspan=1)  # axes_button_left
ax_br = plt.subplot2grid((5,5),(4,3),colspan=2,rowspan=1)  # axes_button_right

# create forward/backward buttons
butt_l = Button(ax_bl, '\N{leftwards arrow}') # or u'' on python 2
butt_r = Button(ax_br, '\N{rightwards arrow}') # or u'' on python 2

# create initial plot
# store index of data and handle to plot as axes property because why not
ax_pl.idat = 0
hplot = ax_pl.scatter(*dat[ax_pl.idat].T)
ax_pl.hpl = hplot
ax_pl.axis('scaled')
ax_pl.axis([dat[...,0].min(),dat[...,0].max(),
            dat[...,1].min(),dat[...,1].max()])
ax_pl.set_autoscale_on(False)
ax_pl.set_title('{}/{}'.format(ax_pl.idat,dat.shape[0]-1))

# define and hook callback for buttons
def replot_data(ax_pl,dat):
    '''replot data after button push, assumes constant data shape'''
    ax_pl.hpl.set_offsets(dat[ax_pl.idat])
    ax_pl.set_title('{}/{}'.format(ax_pl.idat,dat.shape[0]-1))
    ax_pl.get_figure().canvas.draw()

def left_onclicked(event,ax=ax_pl,dat=dat):
    '''try to decrement data index, replot if success'''
    if ax.idat > 0:
        ax.idat -= 1
        replot_data(ax,dat)

def right_onclicked(event,ax=ax_pl,dat=dat):
    '''try to increment data index, replot if success'''
    if ax.idat < dat.shape[0]-1:
        ax.idat += 1
        replot_data(ax,dat)

butt_l.on_clicked(left_onclicked)
butt_r.on_clicked(right_onclicked)

plt.show()