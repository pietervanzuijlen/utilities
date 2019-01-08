'''
Creates environment for localizing the error to elements.

targetdomain: Type of domain on which the values are evaluated
val: namespace based value to be evalueted per element
inttype: type of integration, can be internal, interface or boundary

'''

from nutils import *


class functionbased:

    def __init__(self, domain, geom, basis, degree):
        self.domain = domain
        self.geom   = geom
        self.degree = degree
        self.basis  = basis
        self.type   = 'functionbased'

        self.indicators = {} 

        for i in range(len(basis)):
            self.indicators[i] = 0

    def add(self, targetdomain, basis, val):

        # calculate residual evaluated over each basisfunction
        res = targetdomain.integrate(basis*val*function.J(self.geom), degree=self.degree*2)
        for i in range(len(basis)):
            self.indicators[i] += res[i]

        return self

    def abs_indicators(self):
        return {key:abs(val) for key, val in self.indicaters.items()}
   
class elementbased:

    def __init__(self, domain, geom, degree, dualspacetype=None):
        self.domain = domain
        self.geom   = geom
        self.degree = degree
        self.indicators = {} 
        self.elemsizes  = {}
        self.type   = 'elementbased'
        self.dualspacetype  = dualspacetype

        # get element sizes
        sizes = self.domain.integrate_elementwise(function.J(self.geom), degree=self.degree*2)

        # initiate indicators dictionary 
        for elem, h_K in zip(self.domain, sizes):
            head, tail = transform.lookup_item(elem.transform, self.domain.edict)
            self.indicators[head] = 0
            self.elemsizes[head]  = h_K

        assert len(self.elemsizes) == len(self.domain), 'Amount of element areas is not equal to amount of elements'
        assert len(self.indicators) == len(self.domain), 'Amount of indicators is not equal to amount of elements'

 
    def residualbased(self, targetdomain, val, inttype):
    
        # calculate the residual contribution per element
        res   = targetdomain.integrate_elementwise(val*function.J(self.geom), degree=self.degree*3)
        
        # assert variable lengths
        assert len(res) == len(targetdomain), 'Length of residual and target domain are not equal'
        assert inttype in ['internal','interface','boundary'], 'unvalid inttype'
    
        # factor depends on dimension of integration
        if inttype == 'internal':
            a = 1 
        else:
            a = .5
    
        # link residuals to elements
        for elem, ires in zip(targetdomain, res):
            head, tail = transform.lookup_item(elem.transform, self.domain.edict)
            self.indicators[head] += ires * self.elemsizes[head]**a
            if inttype == 'interface':
                head, tail = transform.lookup_item(elem.opposite, self.domain.edict)
                self.indicators[head] += ires * self.elemsizes[head]**a
    
        return self

    
    def goaloriented(self, targetdomain, val, inttype):

        assert self.dualspacetype == 'p-refined' or self.dualspacetype == 'k-refined', 'The dual space should be enriched with either `p-refined` or `k-refined`. Specify the enrichment method'
    
        # calculate the residual contribution per element
        res   = targetdomain.integrate_elementwise(val*function.J(self.geom), degree=self.degree*2)
        
        # assert variable lengths
        assert len(res) == len(targetdomain), 'Length of residual and target domain are not equal'
        assert inttype in ['internal','interface','boundary'], 'unvalid inttype'
    
        # link residuals to elements
        for elem, ires in zip(targetdomain, res):
            if self.dualspacetype == 'p-refined':
                head, tail = transform.lookup_item(elem.transform[:-1], self.domain.edict)
            elif self.dualspacetype == 'k-refined':
                head, tail = transform.lookup_item(elem.transform, self.domain.edict)
            self.indicators[head] += ires
            if inttype == 'interface':
                if self.dualspacetype == 'p-refined':
                    head, tail = transform.lookup_item(elem.opposite[:-1], self.domain.edict)
                elif self.dualspacetype == 'k-refined':
                    head, tail = transform.lookup_item(elem.opposite, self.domain.edict)
                self.indicators[head] += ires

        return self

    
    def abs_indicators(self):
        return {key:abs(val) for key, val in self.indicators.items()}
