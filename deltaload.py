#! /usr/bin/python3

from nutils import *
import numpy, json, matplotlib.collections

def support_by_dof(domain, basis):
  index = function.concatenate([indices[0] for indices, func in function.blocks(basis)])
  supp = [[] for i in range(len(basis))]
  for trans in domain.transforms:
    for i in index.eval(_transforms=[trans])[0]:
      supp[i].append(trans)
  return supp

def main(poissonsratio = .499,
         youngsmodulus = 3.,
         elasticity = 'neo',
         revolved = False,
         basistype = 'th-std',
         degree = 2,
         radius = 175.,
         thickness = 50.,
         surftens = 35.,
         surftens_inner = 35.,
         surftens_outer = 35.,
         pressure_inner = 0.,
         pressure_outer = 0.,
         refine = .1):

  contactangle = numpy.arccos((surftens_outer-surftens_inner) / surftens)
  neumannangle = numpy.arccos((surftens**2-surftens_inner**2-surftens_outer**2) / (2*surftens_inner*surftens_outer))
  log.user('contact angle: {:.2f} deg'.format(contactangle*180/numpy.pi))
  log.user('neumann angle: {:.2f} deg'.format(neumannangle*180/numpy.pi))

  domain2d, geom2d = mesh.rectilinear([70,10])
  geom2d *= thickness/10
  angles = {}, {}
  istep = 0
  fine = False

  while True:
   with log.context('fine' if fine else 'coarse', istep):

    domain = domain.refined if fine else domain2d * topology.RevolutionTopology() if revolved else domain2d
  
    ns = function.Namespace()
    ns.π = numpy.pi
    ns.dbasis = domain.basis(basistype, degree=degree).vector(2)
    ns.dr = 'dbasis_n0 ?dofs_n'
    ns.dy = 'dbasis_n1 ?dofs_n'
    ns.r, ns.y = ns.x = function.bifurcate1(geom2d) if revolved else geom2d
    if revolved:
      ns.θ = function.bifurcate2(function.RevolutionAngle())
      ns.x_i = '<r cos(θ), y, r sin(θ)>_i'
      ns.dbasis_ni = '<dbasis_n0 cos(θ), dbasis_n1, dbasis_n0 sin(θ)>_i'
    ns.φ = contactangle
    ns.λ = (poissonsratio*youngsmodulus) / ((1+poissonsratio)*(1-2*poissonsratio))
    ns.μ = youngsmodulus / (2*(1+poissonsratio))
    ns.d_i = 'dbasis_ni ?dofs_n'
    ns.X_i = 'x_i + d_i'
    ns.N = ns.X.normal()
    if elasticity == 'neo':
      ns.F_ij = 'X_i,j'
      ns.I = 'F_ij F_ij'
      ns.J = function.determinant(ns.F)
      ns.E = 'λ log(J)^2 / 2 - μ log(J) + μ (I - 3) / 2'
    elif elasticity == 'svk':
      ns.ε_ij = '.5 (d_i,j + d_j,i + d_k,i d_k,j)'
      ns.E = '.5 λ ε_ii ε_jj + μ ε_ij ε_ij'
    elif elasticity == 'lin':
      ns.ε_ij = '.5 (d_i,j + d_j,i)'
      ns.E = '.5 λ ε_ii ε_jj + μ ε_ij ε_ij'
    else:
      raise Exception('invalid elasticity: {!r}'.format(elasticity))
    inner, outer = function.partition(ns.r, radius)
    ns.P = pressure_inner * inner + pressure_outer * outer
    ns.tanhfuncsurften = surftens_inner * inner + surftens_outer * outer
    ns.surftens = surftens
    ns.dXdr_i = 'X_i,0'
    ns.slope = function.arctan2(ns.dXdr[1], ns.dXdr[0])

    sqr = domain.boundary['right,left,bottom'].integral('dr^2' @ ns, degree=degree*2) + domain.boundary['bottom'].integral('dy^2' @ ns, degree=degree*2)
    cons = solver.optimize('dofs', sqr, droptol=1e-15)

    ifaces = domain.boundary['top'].interfaces.sample('gauss', 1)
    ring = ifaces.subset(numpy.equal(ifaces.eval(ns.x[0]), radius))
    assert ring.npoints == 1

    energy = domain.integral('E d:x' @ ns, degree=8)
    energy += domain.boundary['top'].integral('tanhfuncsurften d:X' @ ns, degree=8)
    energy += ring.integral('surftens (cos(φ) dr - sin(φ) dy)' @ ns)
    if pressure_inner == pressure_outer:
      energy -= domain.integral('P d:X' @ ns, degree=8)
      res = energy.derivative('dofs')
    else:
      res = energy.derivative('dofs')
      res -= domain.boundary['top'].integral('dbasis_ni P N_i d:X' @ ns, degree=8)
    goal = ring.integral('π + [slope]' @ ns)
  
    if not fine:

      d0 = nsc.d if istep else 0
      sqr = domain.integral(((ns.d - d0)**2).sum(0), degree=degree*2)
      lhs = solver.optimize('dofs', sqr, constrain=cons) # initial condition

      try:
        for lhs, info in log.iter('iter', solver.newton('dofs', res, constrain=cons, lhs0=lhs)):
          log.info('residual: {:.2e}'.format(info.resnorm))
          if info.resnorm < 1e-8:
            break
      except Exception as e:
        log.error('{} at residual {:.2e}'.format(e, info.resnorm))

      angle = goal.eval(dofs=lhs)
      title = 'angle: {:.2f} deg'.format(angle * 180 / numpy.pi)
      log.user(title)
      angles[fine][len(lhs)] = angle

      bezier = domain.sample('bezier', 5)
      hull = domain.topo1.sample('bezier', 5).hull if revolved else bezier.hull
      X = bezier.eval(ns.X[:2], dofs=lhs)
      bezier = domain.boundary['top'].sample('bezier', 5)
      tri = domain.boundary['top'].topo1.sample('bezier', 5).tri if revolved else bezier.tri
      rsurf, Xsurf, slope = bezier.eval(['x_0', 'X_i', 'slope'] @ ns, dofs=lhs)
      xc, = ring.eval(ns.X, dofs=lhs)

      export.triplot('deformed.jpg', X, hull=hull)

      with export.mplfigure('closeup.jpg') as fig:
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.add_collection(matplotlib.collections.LineCollection(X[hull], linewidths=.1, colors='k'))
        ax.add_collection(matplotlib.collections.LineCollection(Xsurf[tri,:2], colors='k'))
        ax.set_xlim(radius-40, radius+40)
        ax.set_ylim(0, 60)
    
      with export.mplfigure('supercloseup.jpg') as fig:
        ax = fig.add_subplot(111, xlim=(xc[0]-.5, xc[0]+.5), ylim=(xc[1]-.6, xc[1]+.2))
        ax.add_collection(matplotlib.collections.LineCollection(X[hull], linewidths=.1, colors='k'))
        ax.add_collection(matplotlib.collections.LineCollection(Xsurf[tri,:2], colors='k'))
        ax.arrow(xc[0], xc[1], -.1*numpy.cos(contactangle), .1*numpy.sin(contactangle), length_includes_head=True, width=.005, fc='k', ec='k')
  
      with export.mplfigure('angle.jpg') as fig:
        ax = fig.add_subplot(111, title=title)
        ax.add_collection(matplotlib.collections.LineCollection(numpy.array([rsurf, slope*180/numpy.pi]).T[tri]))
        ax.axvline(x=radius, ls=':', color='k')
        ax.axhline(y=0, ls=':', color='k')
        ax.autoscale(enable=True, tight=True)
  
      with export.mplfigure('convergence.jpg') as fig:
        ax = fig.add_subplot(111, xlabel='#dofs', ylabel='neumann angle error [deg]')
        for i in range(2):
          if angles[i]:
            d_, a_ = numpy.array(sorted(angles[i].items())).T
            ax.loglog(d_, abs(a_-neumannangle) * 180 / numpy.pi, 'o*'[i] + ':', label='fine' if i else 'coarse')
        ax.legend()

      z = res.derivative('dofs').eval(dofs=lhs).T.solve(goal.derivative('dofs').eval(dofs=lhs), constrain=numpy.isfinite(cons))
      nsc = ns(dofs=lhs)
      fine = True

    else:

      FF, FC = domain.integrate([(ns.dbasis[:,numpy.newaxis] * ns.dbasis).sum(-1), (ns.dbasis[:,numpy.newaxis] * nsc.dbasis).sum(-1)], degree=degree*2)
      lhsc = FF.solve(FC.matvec(lhs))
      zc = FF.solve(FC.matvec(z))
      z = res.derivative('dofs').eval(dofs=lhsc).T.solve(goal.derivative('dofs').eval(dofs=lhsc), constrain=numpy.isfinite(cons)) #- zc

      bezier = domain.sample('bezier', 5)
      tri = domain.topo1.sample('bezier', 5).tri if revolved else bezier.tri
      x, dual = bezier.eval([ns.x[:2], ns.dbasis.dot(z)])
      with export.mplfigure('dual.png') as fig:
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        im = ax.tripcolor(x[:,0], x[:,1], tri, numpy.linalg.norm(dual, axis=1), shading='gouraud', cmap='jet', rasterized=True, norm=matplotlib.colors.LogNorm())
        fig.colorbar(im, orientation='horizontal')
        ax.autoscale(enable=True, axis='both', tight=True)

      v = res.eval(dofs=lhsc) # v.dot(zc) == 0
      errorest = v.dot(z)
      log.user('error estimate: {:.1e} ({:.0f}% accurate)'.format(v.dot(z) * 180 / numpy.pi, 100*errorest/(angle-neumannangle)))

      supp = support_by_dof(domain, ns.dbasis)
      myrefine = set()
      vz = numpy.abs(v * z)
      toterror = vz.sum()
      referror = 0
      for idof in numpy.argsort(vz)[::-1]:
        transforms = supp[idof]
        if revolved:
          transforms = [trans.trans1[:-1] for trans in transforms]
        minlen = min(len(trans) for trans in transforms)
        myrefine.update(trans for trans in transforms if len(trans) == minlen)
        referror += vz[idof]
        if referror > refine * toterror:
          break

      domain2d = domain2d.refined_by(myrefine)

      fine = False
      istep += 1

if __name__ == '__main__':
  cli.run( main )
