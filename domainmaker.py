from nutils import *
import numpy as np

def rotategrid(grid, geom, angle):

    # construct rotation matrix
    RotMat = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])

    # get function s[ace
    funcsp = grid.basis('th-spline', degree=1 )

    # get control points
    paramcps = grid.project( geom, onto=funcsp.vector(2), geometry=geom, ischeme='gauss4' ).reshape(2,-1).T

    # apply rotation to the control points
    rotparamcps = RotMat.dot(paramcps.T) 

    # construct the new geometery
    geom = funcsp.vector(2).dot(rotparamcps.ravel()) 

    return grid, geom

def lshape(uref=0, width=1, height=1):

    # define patches by verts           # | patchnumber
    patches = [[0,1,3,4],               # | 0 
               [1,2,4,5]]               # | 1

    # define vert positions             # | vertnumber
    verts   = [[0,-.5*height],          # | 0
               [0,0],                   # | 1
               [-.5*width,0],           # | 2
               [.5*width,-.5*height],   # | 3
               [.5*width,.5*height],    # | 4
               [-.5*width,.5*height]]   # | 5

    # make domain
    domain, geom = mesh.multipatch(patches = patches, patchverts=verts, nelems = {None: 1})

    # assign boundaries
    domain = domain.withboundary(
          left    = 'patch1-top',
          right   = 'patch0-right',
          top     = 'patch1-right',
          bottom  = 'patch0-bottom',
          inner   = 'patch0-left,patch1-left',)

    # apply uniform refinement
    domain = domain.refine(uref)

    return domain, geom

def porous(*args, uref=0, rm=.2, r1=.2, r2=.2, r3=.2, r4=.2, M0=.5, M1=.5, **kwargs):
    
    # define patches by verts       # | patchnumber
    patches=[[ [6, 7], [1, 4]],     # | 0
             [ [3, 6], [0, 1]],     # | 1
             [ [1, 4], [9, 8]],     # | 2
             [ [2, 1], [5, 9]],     # | 3
             [ [9, 8], [14, 11]],   # | 4
             [ [12, 9], [15, 14]],  # | 5
             [ [14, 11], [6, 7]],   # | 6
             [ [13, 14], [10, 6]]]  # | 7

    # define vert positions         # | vertnumber
    patchverts = [[r1,0],           # | 0  
                  [(1+r1-r2)/2,0],  # | 1
                  [1-r2,0],         # | 2
                  [0,r1],           # | 3
                  [M0,M1-rm],       # | 4
                  [1,r2],           # | 5
                  [0,(1+r1-r3)/2],  # | 6
                  [M0-rm,M1],       # | 7
                  [M0+rm,M1],       # | 8
                  [1,(1+r2-r4)/2],  # | 9 
                  [0,1-r3],         # | 10 
                  [M0,M1+rm],       # | 11
                  [1,1-r4],         # | 12
                  [r3,1],           # | 13
                  [(1+r3-r4)/2,1],  # | 14
                  [1-r4,1]]         # | 15
                                    
    # make multipatch               
    domain, param = mesh.multipatch(patches=patches,patchverts=patchverts, nelems={None: 1})
    #domain, param = mesh.multipatch(patches=patches,patchverts=patchverts, nelems={None: 1, (3,0):2, (6,1):2, (7,4):2, (2,5):2, (1,9):2})

    # define the function basis of degree 2 
    funcsp = domain.basis('th-spline', degree=2 )

    # get the initial control points 
    paramcps = domain.project( param, onto=funcsp.vector(2), geometry=param, ischeme='gauss4' ).reshape(2,-1).T

    # assign control point position
    paramcps[5]  = [M0-rm,M1-rm]
    paramcps[11] = [r1,r1]
    paramcps[17] = [M0+rm,M1-rm]
    paramcps[23] = [1-r2,r2]
    paramcps[29] = [M0+rm,M1+rm]
    paramcps[35] = [1-r4,1-r4]
    paramcps[41] = [M0-rm,M1+rm]
    paramcps[44] = [r3,1-r3]

    # assign control points weight
    cws = np.ones(len(funcsp))
    cws[5]  = 1/2**.5
    cws[11] = 1/2**.5
    cws[17] = 1/2**.5 
    cws[23] = 1/2**.5
    cws[29] = 1/2**.5
    cws[35] = 1/2**.5
    cws[41] = 1/2**.5
    cws[44] = 1/2**.5
   
    # make nurbes with new control points and weights
    weightfunc = funcsp * cws.T 
    nurbsfunc  = weightfunc / weightfunc.sum([0]) 
 
    # apply uniform refinement
    domain = domain.refine(uref)

    # define new geometry
    geom = nurbsfunc.vector(2).dot(paramcps.T.ravel()) 
    
    # assign boundaries
    domain = domain.withboundary(
          left    = 'patch1-left,patch7-right',
          right   = 'patch3-right,patch5-left',
          top     = 'patch5-right,patch7-left,',
          bottom  = 'patch1-right,patch3-left',
          corners = 'patch1-bottom,patch3-bottom,patch5-bottom,patch7-bottom',
          circle  = 'patch0-top,patch2-top,patch4-top,patch6-top',
          )

    return domain, geom

def cylinder(*args, uref=0, r=.2, M0=.5, M1=.5, H=1, W=1, **kwargs):

    sr = np.sqrt(.5)*r
    qr = np.sqrt(2)*r
    
    # define patches by verts       # | patchnumber
    patches=[[ [0, 2], [1, 3]],     # | 0
             [ [1, 3], [7, 5]],     # | 1
             [ [7, 5], [6, 4]],     # | 2
             [ [6, 4], [0, 2]]]     # | 3

    # define vert positions         # | vertnumber
    patchverts = [[0,0],            # | 0  
                  [W,0],            # | 1
                  [M0-sr,M1-sr],    # | 2
                  [M0+sr,M1-sr],    # | 3
                  [M0-sr,M1+sr],    # | 4
                  [M0+sr,M1+sr],    # | 5
                  [0,H],            # | 6
                  [W,H]]            # | 7
                                    
    # make multipatch               
    domain, param = mesh.multipatch(patches=patches,patchverts=patchverts, nelems={None: 1})

    # define the function basis of degree 2 
    funcsp = domain.basis('th-spline', degree=2 )

    # get the initial control points 
    paramcps = domain.project( param, onto=funcsp.vector(2), geometry=param, ischeme='gauss4' ).reshape(2,-1).T

    # assign control point position
    paramcps[5]  = [M0,M1-qr]
    paramcps[11] = [M0+qr,M1]
    paramcps[17] = [M0,M1+qr]
    paramcps[23] = [M0-qr,M1]

    # assign control points weight
    cws = np.ones(len(funcsp))
    cws[5]  = 1/2**.5
    cws[11] = 1/2**.5
    cws[17] = 1/2**.5 
    cws[23] = 1/2**.5
   
    # make nurbes with new control points and weights
    weightfunc = funcsp * cws.T 
    nurbsfunc  = weightfunc / weightfunc.sum([0]) 
 
    # apply uniform refinement
    domain = domain.refine(uref)

    # define new geometry
    geom = nurbsfunc.vector(2).dot(paramcps.T.ravel()) 
    
    # assign boundaries
    domain = domain.withboundary(
          left    = 'patch3-bottom',
          right   = 'patch1-bottom',
          top     = 'patch2-bottom',
          bottom  = 'patch0-bottom',
          circle  = 'patch0-top,patch1-top,patch2-top,patch3-top',
          )

    return domain, geom

def complexchannel(maxrefine, nelem=21, cutval=0.7):
    
    grid, geom = mesh.rectilinear([np.linspace(0,2,nelem*2+1), np.linspace(0,1,nelem+1)])
    x, y = geom

    level1 = function.Sin(y*0.83*np.pi + 0.1*x + 0.21*function.Sin(x*1.11*np.pi))
    level2 = 3*function.Sin(y*0.97*np.pi + 1.1*x + 0.4*x**2)-2
    level3 = 0.5*(2*function.Sin(y*np.pi-0.6) + 0.8*function.sqrt(x) - 2.6)
    level4 = 0.5*(function.Sin(x*np.pi*0.6-0.8) * function.Sin(y*0.83*np.pi  + 0.13*x + 0.26*function.Sin(x*1.11*np.pi)))**8

    levelset = function.max(level1,level2)+ function.max(0,level3) - level4

    domain = grid.trim(levelset - 0.7 , maxrefine=maxrefine)

    return grid, geom, domain


def channels(uref=0, W=2, H=1, h1=0.2, h2=0.4, w1=0.2, w2=0.2, elemsize=0.1):

    assert h1/elemsize == round(h1/elemsize), 'mimimal step size for length is elemsize'
    assert h2/elemsize == round(h2/elemsize), 'mimimal step size for length is elemsize'
    assert w1/elemsize == round(w1/elemsize), 'mimimal step size for length is elemsize'
    assert w2/elemsize == round(w2/elemsize), 'mimimal step size for length is elemsize'

    # define patches by verts           # | patchnumber
    patches = [[0,4,1,5],               # | 0 
               [1,5,2,6],               # | 1
               [2,6,3,7],               # | 2
               [4,8,5,9],               # | 3
               [6,10,7,11],             # | 4
               [8,12,9,13],             # | 5
               [9,13,10,14],            # | 6
               [10,14,11,15]]           # | 7

    # define vert positions             # | vertnumber
    verts   = [[0,0],                   # | 0
               [w1,0],                  # | 1
               [W-w2,0],                # | 2
               [W,0],                   # | 3
               [0,h1],                  # | 4
               [w1,h1],                 # | 5
               [W-w2,h1],               # | 6
               [W,h1],                  # | 7
               [0,H-h2],                # | 8
               [w1,H-h2],               # | 9
               [W-w2,H-h2],             # | 10
               [W,H-h2],                # | 11
               [0,H],                   # | 12
               [w1,H],                  # | 13
               [W-w2,H],                # | 14
               [W,H]]                   # | 15

    nelems  = {(0,1):       int((w1)/elemsize),
               (1,2):       int((W-w1-w2)/elemsize),
               (2,3):       int((w2)/elemsize),
               (4,5):       int((w1)/elemsize),
               (5,6):       int((W-w1-w2)/elemsize),
               (6,7):       int((w2)/elemsize),
               (8,9):       int((w1)/elemsize),
               (9,10):      int((W-w1-w2)/elemsize),
               (10,11):     int((w2)/elemsize),
               (12,13):     int((w1)/elemsize),
               (13,14):     int((W-w1-w2)/elemsize),
               (14,15):     int((w2)/elemsize),
               (0,4):       int((h1)/elemsize),
               (4,8):       int((H-h1-h2)/elemsize),
               (8,12):      int((h2)/elemsize),
               (1,5):       int((h1)/elemsize),
               (5,9):       int((H-h1-h2)/elemsize),
               (9,13):      int((h2)/elemsize),
               (2,6):       int((h1)/elemsize),
               (6,10):      int((H-h1-h2)/elemsize),
               (10,14):     int((h2)/elemsize),
               (3,7):       int((h1)/elemsize),
               (7,11):      int((H-h1-h2)/elemsize),
               (11,15):     int((h2)/elemsize)}

    # make domain
    domain, geom = mesh.multipatch(patches = patches, patchverts=verts, nelems=nelems)

    # apply uniform refinement
    domain = domain.refine(uref)

    # assign boundaries
    domain = domain.withboundary(
          left    = 'patch0-left,patch3-left,patch5-left',
          right   = 'patch2-right,patch4-right,patch7-right',
          top     = 'patch5-top,patch6-top,patch7-top',
          bottom  = 'patch0-bottom,patch1-bottom,patch2-bottom',
          center  = 'patch1-top,patch3-right,patch4-left,patch6-bottom',
          )

    return domain, geom

def annulus(*args, uref=0, Rin=1, Rout=4, **kwargs):

    Rmid = (Rin + Rout)/2

    domain, param = mesh.rectilinear([numpy.linspace(0,1,2),numpy.linspace(0,1,2)])

    # define the function basis of degree 2 
    funcsp = domain.basis('th-spline', degree=2 )

    # get the initial control points 
    paramcps = domain.project( param, onto=funcsp.vector(2), geometry=param, ischeme='gauss4' ).reshape(2,-1).T

    # assign control point position
    paramcps[0] = [Rin,0]
    paramcps[1] = [Rin,Rin]
    paramcps[2] = [0,Rin]
    paramcps[3] = [Rmid,0]
    paramcps[4] = [Rmid,Rmid]
    paramcps[5] = [0,Rmid]
    paramcps[6] = [Rout,0]
    paramcps[7] = [Rout,Rout]
    paramcps[8] = [0,Rout]

    # assign control points weight
    cws = np.ones(len(funcsp))
    cws[1] = 1/2**.5
    cws[4] = 1/2**.5
    cws[7] = 1/2**.5 
   
    # make nurbes with new control points and weights
    weightfunc = funcsp * cws.T 
    nurbsfunc  = weightfunc / weightfunc.sum([0]) 
 
    # apply uniform refinement
    domain = domain.refine(uref)

    # define new geometry
    geom = nurbsfunc.vector(2).dot(paramcps.T.ravel()) 

    return domain, geom


class immersed:

    def __init__(self, domain, grid):
        
        # Get indices of elements of the grid which are inside the trimmed domain
        indices = np.array([i for i, trans in enumerate(grid.transforms) if domain.transforms.contains(trans)])

        # Construct the list of references for the background topology
        gridrefs = grid.references
        backgroundrefs = []
        for i, ref in enumerate(gridrefs):
            if i in indices:
                backgroundrefs.append(ref)
            else:
                backgroundrefs.append(ref.empty)
        
        # Define the background topoloty as a subset of the grid
        background = grid[numpy.sort(numpy.fromiter(map(grid.transforms.index, domain.transforms), dtype=int))]
        skeleton = background.interfaces

        # Get indices of elements on the trimmed boundary
        trimindices = [domain.transforms.index_with_tail(trans)[0] for trans in domain.boundary['trimmed'].transforms]
        # Filter duplications
        trimindices = sorted(set(trimindices))
        # Translate indices to indices in the background
        backgroundindi = [background.transforms.index(trans) for trans in domain.transforms[numpy.array(trimindices)]]

        grefs = []
        # Find references of the interfaces in the cutelements
        for iface, oppo, ref in zip(skeleton.transforms, skeleton.opposites, skeleton.references):
            ifacehead = background.transforms.index_with_tail(iface)[0]
            oppohead = background.opposites.index_with_tail(oppo)[0]
            if ifacehead in backgroundindi or oppohead in backgroundindi:
                grefs.append(ref)
            else:
                grefs.append(ref.empty)

        assert len(skeleton) == len(grefs), 'Lengths don`t comply'

        # Define the ghost topoloty as a subset of the skeleton
        ghost = topology.SubsetTopology(background.interfaces, grefs)

        self.background = background
        self.skeleton = skeleton
        self.ghost = ghost

#def skelghost(grid, trim):
#
#    # Get indices of elements of the grid which are inside the trimmed domain
#    indices = np.array([i for i, trans in enumerate(grid.transforms) if trim.transforms.contains(trans)])
#
#    # Construct the list of references for the background topology
#    gridrefs = grid.references
#    backgroundrefs = []
#    for i, ref in enumerate(gridrefs):
#        if i in indices:
#            backgroundrefs.append(ref)
#        else:
#            backgroundrefs.append(ref.empty)
#    
#    # Define the background topoloty as a subset of the grid
##    background = topology.SubsetTopology(grid, backgroundrefs)
#    background = grid[numpy.sort(numpy.fromiter(map(grid.transforms.index, trim.transforms), dtype=int))]
#    skeleton = background.interfaces
#
#    # Get indices of elements on the trimmed boundary
#    trimindices = [trim.transforms.index_with_tail(trans)[0] for trans in trim.boundary['trimmed'].transforms]
#    # Filter duplications
#    trimindices = sorted(set(trimindices))
#    # Translate indices to indices in the background
#    backgroundindi = [background.transforms.index(trans) for trans in trim.transforms[numpy.array(trimindices)]]
#
#    grefs = []
#    # Find references of the interfaces in the cutelements
#    for iface, oppo, ref in zip(skeleton.transforms, skeleton.opposites, skeleton.references):
#        ifacehead = background.transforms.index_with_tail(iface)[0]
#        oppohead = background.opposites.index_with_tail(oppo)[0]
#        if ifacehead in backgroundindi or oppohead in backgroundindi:
#            grefs.append(ref)
#        else:
#            grefs.append(ref.empty)
#
#    assert len(skeleton) == len(grefs), 'Lengths don`t comply'
#
#    # Define the ghost topoloty as a subset of the skeleton
#    ghost = topology.SubsetTopology(background.interfaces, grefs)
#
#    return background.interfaces, ghost
