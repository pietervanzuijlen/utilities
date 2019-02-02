from   nutils import *
import numpy 
import matplotlib.pyplot as plt
from   matplotlib import collections

def plot_indicators(name, domain, geom, indicators, npoints=5, shape=0, bartitle='', alpha=0, normalize=False):

    # Amount of subplots
    n = len(indicators)

    # Define subplot shape
    if shape == 0:
        nsqrt = numpy.sqrt(n)
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
        if vmax <= max(val):
            vmax = max(val)
        # Addition for non-absolute indicaters
        if vmin >= min(val):
            vmin = min(val)
        colors[key] = {} 
        for i, trans in enumerate(domain.transforms):
            colors[key][trans] = val[i]

    # Get domain shape
    bezier = domain.sample('bezier', npoints)
    x = bezier.eval(geom)

    # Export the figure
    with export.mplfigure(name+'.png') as fig:

        # Making the subplots
        i = 0 
        for key, val in indicators.items():
            trans = domain.transforms
            colors = tuple(val) 
            elemval = function.elemwise(trans, colors)
            elemvalue = bezier.eval(elemval)
            i += 1
            subshape = int(shape*10+i)
            ax = fig.add_subplot(subshape, aspect='equal')

            if normalize:
                im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, elemvalue, shading='gouraud', cmap='summer', vmin=vmin, vmax=vmax)
            else:
                im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, elemvalue, shading='gouraud', cmap='summer')
                fig.colorbar(im)
            ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.3, alpha=alpha))
            ax.autoscale(enable=True, axis='both', tight=True)
            ax.set_title(str(key))
            ax.axis('off')

        # Set the colorbar 
        if normalize:
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

    return 

def plot_interfaces(name, domain, geom, interfaces, npoints=5, color=0.5, cmap='summer', title='', alpha=1):

    vmin = 0
    vmax = 1

    domsample = domain.sample('bezier', npoints)
    x_domain, fill = domsample.eval([geom, color])
    ifacesample = interfaces.sample('bezier', npoints*10)
    x_iface = ifacesample.eval(geom)

    with export.mplfigure(name+'.png') as fig:
      ax = fig.add_subplot(111, aspect='equal')
      im = ax.tripcolor(x_domain[:,0], x_domain[:,1], domsample.tri, fill, shading='gouraud', cmap=cmap, alpha=alpha)
      ax.scatter(x_iface[:,0],x_iface[:,1],c='r',s=.8)

    return 

def plot_levels(name, domain, geom, minlvl=1, npoints=5, cmap='summer', title='', alpha=1):

    levels = numpy.array([len(trans) for trans in domain.transforms])
    levels = tuple(levels-minlvl-1)
    lvl = function.elemwise(domain.transforms, levels)

    bezier = domain.sample('bezier', npoints)
    x,fill = bezier.eval([geom,lvl])

    # background color is found as the color value between 0 and 1 projected on the colormap
    with export.mplfigure(name+'.png') as fig:
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, fill, shading='gouraud', cmap=cmap, alpha=alpha)
        ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5))
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)

        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title='Levels')
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')

        ax.axis('off')

    return 

def plot_solution(name, domain, geom, val, npoints=5, cmap='jet', title='', bartitle='', alpha=0, grid=False, **arguments):

    bezier = domain.sample('bezier', npoints)
    x,val = bezier.eval([geom,val], arguments=arguments)

    with export.mplfigure(name+'.png') as fig:
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, val, shading='gouraud', cmap=cmap)
        # Add grid mesh instead of regular mesh
        if grid:
            bezier = grid.sample('bezier', npoints)
            x = bezier.eval(geom)
        ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5, alpha=alpha))
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)
        
        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title=bartitle)
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')

    return


def plot_streamlines(name, domain, geom, ns, val, npoints=5, cmap='jet', title='', bartitle='', every=.05, alpha=0, grid=False, **arguments):

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
        ax.tricontour(x[:,0], x[:,1], bezier.tri, stream, 16, colors='k', linestyles='dotted', linewidths=.5, zorder=9)
        # Add grid mesh instead of regular mesh
        if grid:
            bezier = grid.sample('bezier', npoints)
            x = bezier.eval(geom)
        ax.add_collection(collections.LineCollection(x[bezier.hull], colors='w', linewidths=.5, alpha=alpha))
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)

        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title=bartitle)
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')
  
    return


def plot_convergence(name, xval, yval, labels=None, title='', levels={}, slopemarker=None):
    
    from mpltools import annotation

    assert yval.keys() == xval.keys(), 'the axes keys should have the same names'
   
    color = ['c','b','y','g','y']

    levelxpts = {}
    levelypts = {}

    with export.mplfigure(name+'.png') as fig:
        ax = fig.add_axes([.2,.1,.75,.8])
        #ax = fig.add_subplot(111)
        ax.autoscale(enable=True, axis='both', tight=True)

        # Add labels if given
        if labels:
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])

        for i, key in enumerate(yval.keys()):
    
            # make convergence plot
            im = ax.loglog(xval[key], yval[key], color[i])

            # Indicating hierarchical levels
            if key in levels:
                vals = levels[key]
                levelxpts[key] = []
                levelypts[key] = []
                lvl = vals[0] 
                for j, val in enumerate(vals):
                    if val > lvl:
                        levelxpts[key] += [xval[key][j]]
                        levelypts[key] += [yval[key][j]]
                        lvl = val 
                im = ax.scatter(levelxpts[key], levelypts[key], marker="+", s=8**2, c=color[i])

            # Add slope indicator if given
            if slopemarker:
                
               # xlog = numpy.log10(xval[key])
               # ylog = numpy.log10(yval[key])

               # fit  = numpy.polyfit(xlog,ylog,1)
               # log.user(fit)
               # 
               # slope = slopemarker[key][0]
               # dist  = slopemarker[key][1]

               # xmin = xval[key][0]
               # xmax = xval[key][-1]
               # xpos = 10**((numpy.log10(xmax)+numpy.log10(xmin))/2)

               # ymin = yval[key][0]
               # ymax = yval[key][-1]
               # ypos = 10**((numpy.log10(ymax)+numpy.log10(ymin))/(2-dist))

               # poly_settings = {'edgecolor': color[i]}

               # annotation.slope_marker((xpos,ypos), slope, ax=ax, poly_kwargs=poly_settings, invert=True)

                L = 0.3
                
                xlog = numpy.log10(xval[key])
                ylog = numpy.log10(yval[key])

                fit  = numpy.polyfit(xlog,ylog,1)
                log.user(fit)

                slope = int(round(fit[0]))
                
                xmin = xval[key][0]
                xmax = xval[key][-1]
                Lx   = numpy.log10(xmax)-numpy.log10(xmin)

                x1   = 10**(numpy.log10(xmin)+Lx*(0.5-L/2))
                x2   = 10**(numpy.log10(xmin)+Lx*(0.5+L/2))
                xtxt = 10**(numpy.log10(xmin)+Lx*0.5)

                ymin = yval[key][0]
                ymax = yval[key][-1]
                Ly   = numpy.log10(ymax)-numpy.log10(ymin)

                y1   = 10**(numpy.log10(ymin)+Ly*(0.6-L/2))
                y2   = 10**(numpy.log10(ymin)+Ly*(0.6-L/2)+Lx*slope*L)
                ytxt = 10**(numpy.log10(ymin)+Ly*(0.6-L/2)+Lx*slope*L/2)
                
                im = ax.loglog([x1,x2], [y1,y2] , 'k')
                im = ax.text(xtxt,ytxt,str(slope))

        # Add legend and title
        ax.legend(yval.keys())
        ax.set_title(title)

   
    return



