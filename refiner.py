'''
uses the indicators defined with the module indicater to refine the domain. The type of indicator is automatically determined with the indicator.type value. The domain can be refined with either a fixed amount of to be refined objects or with a fraction of the error value.

domain: domain to be refined
indicators: object constructed with the indicater module
num: defines how much to refine 
* if num is a int, that amount of objects will be refined
* if num is a float between 0 and 1, all objects which have an indicator higher than (1-num)*max(indicator) will be refined 

'''

from nutils import *
import numpy

def refine(domain, indicators, num, maxlevel=10):

    if type(num) == int:
        assert len(domain) >= num, 'Amount of functions to be refined should be lower than the total amount of elements'
        threshold = sorted(list(indicators.indicators.values()))[-num]
    elif type(num) == float:
        assert num < 1, 'Fraction should be lower than 1'
        threshold = max(list(indicators.indicators.values()))*(1-num)

    to_refine = []
    
    for key in list(indicators.indicators.keys()):
        if indicators.indicators[key] >= threshold:
            to_refine += [key]
            
    if indicators.type == 'elementbased':
        assert len(domain) == len(indicators.indicators), 'Amount of indicators and amount of elements are not equal'
        marked = numpy.array(domain.elements)[to_refine]
    elif indicators.type == 'functionbased':
        assert len(indicators.basis) == len(indicators.indicators), 'Amount of indicators and amount of basis functions are not equal'
        mask = numpy.zeros(len(indicators.basis), dtype=bool)
        mask[to_refine] = True
        marked = domain.supp(indicators.basis, mask)

    for elem in marked:
        if len(elem.transform) <= maxlevel+1:
            domain = domain.refined_by((elem.transform,))

    return domain
