from   nutils import *
import numpy 
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_indicators(name, domain, geom, indicators, npoints=5, shape=0, bartitle='', alpha=0, normalize=False, labelsize=10):

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
        fig.patch.set_alpha(0)

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
                cbar_ax = fig.add_axes([0.85, 0.11, 0.03, 0.77], title=bartitle)
                fig.colorbar(im, cax=cbar_ax)
                cbar_ax.tick_params(labelsize=labelsize)
            ax.add_collection(mpl.collections.LineCollection(x[bezier.hull], colors='k', linewidths=.3, alpha=alpha))
            ax.autoscale(enable=True, axis='both', tight=True)
            ax.set_title(str(key))
            ax.axis('off')

        # Set the colorbar 
        if normalize:
            fig.subplots_adjust(left=0.05, right=0.75)
            cbar_ax = fig.add_axes([0.85, 0.11, 0.03, 0.77], title=bartitle)
            fig.colorbar(im, cax=cbar_ax)
            cbar_ax.tick_params(labelsize=labelsize)

    return 

def plot_mesh(name, domain, geom, npoints=5, color=0, cmap='Greys', title='', axes=False, figsize=[6,6], linewidths=.5):

    vmin = 0
    vmax = 1

    bezier = domain.sample('bezier', npoints)
    x,fill = bezier.eval([geom,color])

    # background color is found as the color value between 0 and 1 projected on the colormap
    with export.mplfigure(name+'.png') as fig:
        fig.patch.set_alpha(0)
        fig.set_size_inches(figsize[0],figsize[1])
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, fill, shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.add_collection(mpl.collections.LineCollection(x[bezier.hull], colors='k', linewidths=linewidths))
        ax.set_title(title)
        if axes == False:
            ax.axis('off')

    return 

def plot_edge(name, domain, geom, interfaces, npoints=5, color=0.5, cmap='summer', title='', alpha=1):

    vmin = 0
    vmax = 1

    domsample = domain.sample('bezier', npoints)
    x_domain, fill = domsample.eval([geom, color])
    ifacesample = interfaces.sample('bezier', npoints*10)
    x_iface = ifacesample.eval(geom)

    with export.mplfigure(name+'.png') as fig:
      fig.patch.set_alpha(0)
      ax = fig.add_subplot(111, aspect='equal')
      im = ax.tripcolor(x_domain[:,0], x_domain[:,1], domsample.tri, fill, shading='gouraud', cmap=cmap, alpha=alpha)
      ax.scatter(x_iface[:,0],x_iface[:,1],c='r',s=.8)

    return 

def plot_trim(name, domain, background, geom, npoints=5, color=0.5, bcolor=0.2, cmap='viridis', title=''):
    bezierdom = domain.sample('bezier', npoints)
    bezierbg = background.sample('bezier', npoints)
    xdom = bezierdom.eval(geom)
    xbg = bezierbg.eval(geom)
    with export.mplfigure(name+'.png', dpi=150) as fig:
        fig.patch.set_alpha(0)
        fig.set_size_inches(6,6)
        ax = fig.add_axes([0,0,1,1], aspect='equal')
        ax.add_collection(mpl.collections.PolyCollection(xbg[bezierbg.tri], edgecolors='none', facecolors='#ccccccff', antialiaseds=False))
        ax.add_collection(mpl.collections.LineCollection(xbg[bezierbg.hull], colors='w', linewidths=1, antialiaseds=True))
        ax.add_collection(mpl.collections.PolyCollection(xdom[bezierdom.tri], edgecolors='none', facecolors='#3e7ca8ff', antialiaseds=False))
        ax.add_collection(mpl.collections.LineCollection(xdom[bezierdom.hull], colors='w', linewidths=1, antialiaseds=True))
        ax.autoscale_view()
        ax.set_title(title)
        ax.axis('off')

def plot_levels(name, domain, geom, uref=0, npoints=5, cmap='summer', title='', alpha=1):

    levels = numpy.array([len(trans) for trans in domain.transforms])
    maxlevel = max(levels)-uref-1
    levels = tuple(levels-uref-1)
    lvl = function.elemwise(domain.transforms, levels)

    bezier = domain.sample('bezier', npoints)
    x,fill = bezier.eval([geom,lvl])

    colormap = plt.get_cmap(cmap, maxlevel)

    # background color is found as the color value between 0 and 1 projected on the colormap
    with export.mplfigure(name+'.png') as fig:
        fig.patch.set_alpha(0)
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, fill, shading='gouraud', cmap=colormap, alpha=alpha)
        ax.add_collection(mpl.collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5))
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)

        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title='Levels')
        fig.colorbar(im, cax=cbar_ax, ticks=numpy.arange(maxlevel)+1)
        cbar_ax.yaxis.set_ticks_position('right')

        ax.axis('off')

    return 

def plot_solution(name, domain, geom, val, npoints=5, cmap='viridis', title='', bartitle='', alpha=0, axisoff=False, grid=False, vmax=None, vmin=None, figsize=[5,4], cbarsize=[0.8, 0.1, 0.03, 0.8], axsize=[0.01,0.1,.8,.8], labelsize=10, **arguments):

    bezier = domain.sample('bezier', npoints)
    x,val = bezier.eval([geom,val], arguments=arguments)

    with export.mplfigure(name+'.png') as fig:
        fig.patch.set_alpha(0)
        fig.set_size_inches(figsize[0],figsize[1])
        #ax = fig.add_axes([0,0.1,.8,.8], aspect='equal')
        ax = fig.add_axes(axsize, aspect='equal')
        ax.set_xmargin(0)
        ax.set_ymargin(0)
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, val, shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax)
        # Add grid mesh instead of regular mesh
        if grid:
            bezier = grid.sample('bezier', npoints)
            x = bezier.eval(geom)
        ax.add_collection(mpl.collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5, alpha=alpha))
        ax.set_title(title)
        
        #cbar_ax = fig.add_axes([0.8, 0.1, 0.03, 0.8], title=bartitle)
        cbar_ax = fig.add_axes(cbarsize, title=bartitle)
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')
        cbar_ax.tick_params(labelsize=labelsize)

        if axisoff: 
            ax.axis('off')

    return


def plot_streamlines(name, domain, geom, ns, val, npoints=5, cmap='viridis', title='', bartitle='', every=.05, alpha=0, axisoff=False, grid=False, vmax=None, vmin=None, figsize=[5,4], cbarsize=[0.85, 0.1, 0.03, 0.8], axsize=[0.08,0.1,.8,.8], labelsize=10, **arguments):

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
        fig.patch.set_alpha(0)
        fig.set_size_inches(figsize[0],figsize[1])
        ax = fig.add_axes(axsize, aspect='equal')
        #ax = fig.add_axes([0,0.1,.8,.8], aspect='equal')
        ax.set_xmargin(0)
        ax.set_ymargin(0)
        im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, vector, shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.tricontour(x[:,0], x[:,1], bezier.tri, stream, 16, colors='k', linestyles='dotted', linewidths=.5, zorder=9)
        # Add grid mesh instead of regular mesh
        if grid:
            bezier = grid.sample('bezier', npoints)
            x = bezier.eval(geom)
        ax.add_collection(mpl.collections.LineCollection(x[bezier.hull], colors='w', linewidths=.5, alpha=alpha))
        ax.set_title(title)

        #cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title=bartitle)
        cbar_ax = fig.add_axes(cbarsize, title=bartitle)
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')
        cbar_ax.tick_params(labelsize=labelsize)

        if axisoff: 
            ax.axis('off')
  
    return

def plot_contour(name, domain, geom, val, npoints=5, zorder=9, cmap='viridis', title='', bartitle='', alpha=1, tol=1e-10, linewidths=.5, **arguments):

    bezier = domain.sample('bezier', npoints)
    x,val = bezier.eval([geom,val], arguments=arguments)

    val[abs(val) < tol] = 0.0

    with export.mplfigure(name+'.png') as fig:
        fig.patch.set_alpha(0)
        ax = fig.add_axes([.1,.1,.8,.8], aspect='equal')
        im = ax.tricontour(x[:,0], x[:,1], val, 10, linewidths=linewidths, zorder=zorder, cmap=cmap)
        ax.add_collection(mpl.collections.LineCollection(x[bezier.hull], colors='k', linewidths=1, alpha=alpha ))
        ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_ymargin(0)
        ax.axis('off')
        
        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8], title=bartitle)
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar_ax.yaxis.set_ticks_position('right')

    return

def plot_convergence(name, xval, yval, labels=None, marker='o', title='', levels={}, slopemarker=None, fontsize='large'):

    # fontsize 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
    
    #from mpltools import annotation

    assert yval.keys() == xval.keys(), 'the axes keys should have the same names'
   
    #color = ['c','b','y','g','y']
    color = ['orange','crimson','indigo','orchid','tomato']

    levelxpts = {}
    levelypts = {}
    handles = []

    with export.mplfigure(name+'.png') as fig:
        fig.patch.set_alpha(0)
        ax = fig.add_axes([.2,.1,.75,.8])
        #ax = fig.add_subplot(111)
        ax.autoscale(enable=True, axis='both', tight=True)

        # Add labels if given
        if labels:
            ax.set_xlabel(labels[0], fontsize=fontsize)
            ax.set_ylabel(labels[1], fontsize=fontsize)

        for i, key in enumerate(yval.keys()):
    
            # make convergence plot
            im = ax.loglog(xval[key], yval[key], color[i], marker=marker)
            handles += [mpl.lines.Line2D([], [],color=color[i], label=key)]
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
        ax.legend(handles=handles, fontsize=fontsize)
        ax.set_title(title)

   
    return



