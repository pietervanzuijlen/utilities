from   nutils import *
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import collections

def plot_indicators(name, domain, geom, indicators, npoints=5, shape=0, bartitle=''):

    # Amount of subplots
    n = len(indicators)

    # Define subplot shape
    if shape == 0:
        nsqrt = np.sqrt(n)
        if nsqrt <= round(nsqrt):
            shape = round(nsqrt)*10+round(nsqrt)
        else:
            shape = (round(nsqrt)+1)*10+(round(nsqrt)+1)
    else:
        assert n <= int(str(shape)[0])*int(str(shape)[1]), 'Given shape is not large enough'

    # Minimum and maximum value for plotting
    vmin = 0
    vmax = 0 
    colors = {}

    # Assemble colors for the plots
    for key, val in indicators.items():
        if vmax <= max(val.values()):
            vmax = max(val.values())
        # Addition for non-absolute indicaters
        if vmin >= min(val.values()):
            vmin = min(val.values())
        colors[key] = np.array([]) 
        for i in range(len(domain)):
            colors[key] = np.append(colors[key],np.ones(npoints**2)*val[i])

    # Get domain shape
    bezier = domain.sample('bezier', npoints)
    x = bezier.eval(geom)

    # Export the figure
    with export.mplfigure(name+'.png') as fig:

        # Making the subplots
        i = 0 
        for key, val in colors.items():
            i += 1
            subshape = int(shape*10+i)
            ax = fig.add_subplot(subshape, aspect='equal')
            im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, val, shading='gouraud', cmap='summer', vmin=vmin, vmax=vmax)
            ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.3))
            ax.autoscale(enable=True, axis='both', tight=True)
            ax.set_title(str(key))
            ax.axis('off')

        # Set the colorbar 
        fig.subplots_adjust(left=0.05, right=0.75)
        cbar_ax = fig.add_axes([0.85, 0.11, 0.03, 0.77], title=bartitle)
        fig.colorbar(im, cax=cbar_ax)

    return 


def plot_mesh(name, domain, geom, npoints=5, color=0.5, cmap='jet', title=''):

    vmin = 0
    vmax = 1

    bezier = domain.sample('bezier', npoints)
    x,fill = bezier.eval([geom,color])

    # background color is found as the color value between 0 and 1 projected on the colormap
    with export.mplfigure(name+'.png') as fig:
        fig.set_size_inches(6,6)
        ax = fig.add_axes([0,0,1,1], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, fill, shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5))
        ax.set_title(title)

        ax.axis('off')


def plot_solution(name, domain, geom, val, npoints=5, cmap='jet', title='', bartitle='', alpha=0):

    bezier = domain.sample('bezier', npoints)
    x,val = bezier.eval([geom,val])

    with export.mplfigure(name+'.png') as fig:
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, val, shading='gouraud', cmap=cmap)
        ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5, alpha=alpha))
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)
        
        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title=bartitle)
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')

    return


def plot_streamlines(name, domain, geom, ns, val, npoints=5, cmap='jet', title='', bartitle='', every=.05, alpha=0, **arguments):

    ns.vector = val

    # compute streamlines
    ns = ns.copy_()
    ns.streambasis = domain.basis('th-spline', degree=2)[1:] # remove first dof to obtain non-singular system
    ns.stream = 'streambasis_n ?streamdofs_n'
    sqr = domain.integral('(vector_0 - stream_,1)^2 + (vector_1 + stream_,0)^2 d:x' @ ns, degree=4)
    arguments['streamdofs'] = solver.optimize('streamdofs', sqr, arguments=arguments)
  
    # evaluate the values
    bezier = domain.sample('bezier', npoints)
    x, vector, stream = bezier.eval([ns.x, function.norm2(ns.vector), ns.stream], arguments=arguments)

    # plot vector as field and streamlines as dashed
    with export.mplfigure(name+'.png') as fig:
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, vector, shading='gouraud', cmap=cmap)
        ax.add_collection(collections.LineCollection(x[bezier.hull], colors='w', linewidths=.5, alpha=alpha))
        ax.tricontour(x[:,0], x[:,1], bezier.tri, stream, 16, colors='k', linestyles='dotted', linewidths=.5, zorder=9)
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)

        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title=bartitle)
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')
  
    return

