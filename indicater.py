'''
Creates environment for localizing the error to elements.

targetdomain: Type of domain on which the values are evaluated
val: namespace based value to be evalueted per element
inttype: type of integration, can be internal, interface or boundary

'''

from nutils import *



class define:

    def __init__(self, domain, geom, degree):
        self.domain = domain
        self.geom   = geom
        self.degree = degree
        self.indicators = {} 
        self.elemsizes  = {}

        # get element sizes
        sizes = self.domain.integrate_elementwise(function.J(self.geom), degree=self.degree)

        # initiate indicators dictionary 
        for elem, h_K in zip(self.domain, sizes):
            head, tail = transform.lookup_item(elem.transform, self.domain.edict)
            self.indicators[head] = 0
            self.elemsizes[head]  = h_K

        assert len(self.elemsizes) == len(self.domain), 'Amount of element areas is not equal to amount of elements'
        assert len(self.indicators) == len(self.domain), 'Amount of indicators is not equal to amount of elements'

 
    def residualbased(self, targetdomain, val, inttype):
    
        # calculate the residual contribution per element
        res   = targetdomain.integrate_elementwise(val, degree=self.degree)
        
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
            head, tail = transform.lookup_item(elem.transform, self.domain.edict )
            self.indicators[head] += ires * self.elemsizes[head]**a
            if inttype == 'interface':
                head, tail = transform.lookup_item(elem.opposite, self.domain.edict )
                self.indicators[head] += ires * self.elemsizes[head]**a
    
        return self

    
    def goaloriented(self, targetdomain, val, inttype):
    
        # calculate the residual contribution per element
        res   = targetdomain.integrate_elementwise(val, degree=self.degree)
        
        # assert variable lengths
        assert len(res) == len(targetdomain), 'Length of residual and target domain are not equal'
        assert inttype in ['internal','interface','boundary'], 'unvalid inttype'
    
        # link residuals to elements
        for elem, ires in zip(targetdomain, res):
            head, tail = transform.lookup_item(elem.transform[:-1], self.domain.edict )
            self.indicators[head] += ires
            if inttype == 'interface':
                head, tail = transform.lookup_item(elem.opposite[:-1], self.domain.edict )
                self.indicators[head] += ires
    
        return self