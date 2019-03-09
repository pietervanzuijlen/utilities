from nutils import *
import numpy as np

def laplace(domain, geom, ns):

    # Cartesian coordinates
    x, y = geom

    # Polar coordinates
    theta  = np.pi*3/4 - function.ArcTan2(x-y,x+y)
    radius = (x**2 + y**2)**.5

    # Exact solution
    ns.uexact  = radius**(2/3) * function.Sin(2/3*theta) 
    
    return ns


def stokes(domain, geom, ns):

    # Cartesian coordinates
    x, y = geom
    
    # Constants
    alpha = 856399/1572864
    omega = np.pi * 3 / 2
    
    # Polar coordinates
    theta  = np.pi*3/4 - function.ArcTan2(x-y,x+y)
    radius = (x**2 + y**2)**.5

    #ns.alpha = 856399/1572864
    #ns.omega = np.pi * 3 / 2
    #ns.theta  = np.pi*3/4 - function.ArcTan2(x+y,-x+y) 
    #ns.radius = (x**2 + y**2)**.5
    #ns.psi = 'sin((1 + alpha) theta) cos(alpha omega) / (1 + alpha) - cos((1 + alpha) theta) - sin((1 - alpha) theta) cos(alpha omega) / (1 - alpha) + cos((1 - alpha) theta)'
    #ns.d1psi = 'cos((1 + alpha) theta) cos(alpha omega) + (1 + alpha) sin((1 + alpha) theta) - cos((1 - alpha) theta) cos(alpha omega) - (1 - alpha) sin((1 - alpha) theta)'    
    #ns.d3psi = '-(1 + alpha)^2 cos((1 + alpha) theta) cos(alpha omega) - (1 + alpha)^3 sin((1 + alpha) theta) + (1 - alpha)^2 cos((1 - alpha) theta) cos(alpha omega) + (1 - alpha)^3 sin((1 - alpha) theta)'    

    # Angular values    
    psi = (function.Sin((1+alpha)*theta)*function.Cos(alpha*omega)/(1+alpha) 
          -function.Sin((1-alpha)*theta)*function.Cos(alpha*omega)/(1-alpha) 
          - function.Cos((1+alpha)*theta) + function.Cos((1-alpha)*theta))

    d1psi = (function.Cos((1+alpha)*theta)*function.Cos(alpha*omega)
            -function.Cos((1-alpha)*theta)*function.Cos(alpha*omega)
            +(1+alpha)*function.Sin((1+alpha)*theta)-(1-alpha)*function.Sin((1-alpha)*theta))

    d3psi = (-(1+alpha)**2*function.Cos((1+alpha)*theta)*function.Cos(alpha*omega) 
            +(1-alpha)**2*function.Cos((1-alpha)*theta)*function.Cos(alpha*omega)
            -(1+alpha)**3*function.Sin((1+alpha)*theta)+(1-alpha)**3*function.Sin((1-alpha)*theta))


    # Velocity magnitudes
    #ns.uexactx = 'radius^(856399 / 1572864) (-(1 + alpha) cos(theta) psi + sin(theta) d1psi)'
    #ns.uexacty = '-radius^(856399 / 1572864) ((1 + alpha) sin(theta) psi + cos(theta) d1psi)'
    #ns.pexact  = '-radius^(-716465 / 1572864) ((1 + alpha)^2 d1psi + d3psi) / (1 - alpha)'

    ns.uexactx = radius**(1-alpha)*(-(1+alpha)*function.Cos(theta)*psi+function.Sin(theta)*d1psi) 
    ns.uexacty = -radius**(1-alpha)*((1+alpha)*function.Sin(theta)*psi+function.Cos(theta)*d1psi)
    ns.pexact  = -radius**(alpha-1)*((1+alpha)**2*d1psi+d3psi)/(1-alpha)
    
    ns.uexact_i = '<uexactx , uexacty>_i'

    ns.psi = psi
    ns.d1psi = d1psi
    ns.d3psi = d3psi
    ns.theta = theta


    return ns
