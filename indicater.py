from nutils import *
import numpy

def integrate(domain, geom, degree, value, targetdomain, interfaces=None, depth=None):

    # Functionbased
    # Use following line when using functionbased refinement:
    # indicators = targetdomain.integrate(basis*value*function.J(geom), degree=degree*2)

    indicators = numpy.zeros(len(domain))
    res = targetdomain.integrate_elementwise(value*function.J(geom), degree=degree*2)

    heads = [domain.transforms.index_with_tail(trans[:depth])[0] for trans in targetdomain.transforms]
    for i, head in enumerate(heads):
        indicators[head] += res[i]

    if interfaces:
        heads = [domain.transforms.index_with_tail(trans[:depth])[0] for trans in targetdomain.opposites]
        for i, head in enumerate(heads):
            indicators[head] += res[i]

    return indicators
