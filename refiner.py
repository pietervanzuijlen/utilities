'''
uses the indicators defined with the module indicater to refine the domain. The type of indicator is automatically determined with the indicator.type value. The domain can be refined with either a fixed amount of to be refined objects or with a fraction of the error value.

domain: domain to be refined
indicators: object constructed with the indicater module
num: defines how much to refine 
* if num is a int, that amount of objects will be refined
* if num is a float between 0 and 1, all objects which have an indicator higher than (1-num)*max(indicator) will be refined 

For function based, support detectors need to be added
'''

from nutils import *
import numpy

def dorfler_marking(domain, indicators, num, maxlevel=10):

    assert num > 0 and num <= 1, 'the fraction should be between 0 and 1'
    assert indicators.type == 'elementbased' or indicators.type == 'functionbased', 'indicater type should be elementbased or functionbased'
    assert isinstance(num, float), 'num should be a float between 0 and 1' 

    if indicators.type == 'elementbased':

        transforms = []
        values = []
        for elem in domain:
            head, tail = transform.lookup_item(elem.transform, domain.edict)
            transforms.append(elem.transform)
            values.append( abs(indicators.indicators[head]) )

        values = numpy.array(values)

        sortindi = numpy.flip(numpy.argsort(values))
        sortvalues = values[sortindi]
        
        thres  = num * sortvalues.sum()
        cumsum = sortvalues.cumsum()
        select = numpy.zeros(len(cumsum),dtype=bool)
        for i, v in enumerate(cumsum):
            select[i] = True
            if v > thres:
                break
        refindi = sortindi[select]
        reftransforms = tuple(transforms[i] for i in refindi if len(transforms[i]) <= maxlevel+1)

    elif indicators.type == 'functionbased':

        values = []
        for i in indicators.indicators.keys():
            values.append( abs(indicators.indicators[i]) )

        values = numpy.array(values)

        sortindi = numpy.flip(numpy.argsort(values))
        sortvalues = values[sortindi]
        
        thres  = num * sortvalues.sum()
        cumsum = sortvalues.cumsum()
        select = numpy.zeros(len(cumsum),dtype=bool)
        for i, v in enumerate(cumsum):
            select[i] = True
            if v > thres:
                break
        refindi = sortindi[select]

        mask = numpy.zeros(len(indicators.basis), dtype=bool)
        mask[refindi] = True

        reftransforms = tuple(elem.transform for elem in domain.supp(indicators.basis, mask) if len(elem.transform) <= maxlevel+1)

    domain = domain.refined_by(reftransforms)

    return domain

def fractional_marking(domain, indicators, num, maxlevel=10):

    assert indicators.type == 'elementbased' or indicators.type == 'functionbased'
    assert isinstance(num, float), 'num should be a float between 0 and 1' 

    values = [abs(val) for val in indicators.indicators.values()]
    threshold = max(values)*(1-num)
    to_refine = []

    for key in list(indicators.indicators.keys()):
        if abs(indicators.indicators[key]) >= threshold:
            to_refine += [key]

    if indicators.type == 'elementbased':
        marked = numpy.array(domain.elements)[to_refine]
    elif indicators.type == 'functionbased':
        mask = numpy.zeros(len(indicators.basis), dtype=bool)
        mask[to_refine] = True
        marked = domain.supp(indicators.basis, mask)

    for elem in marked:
        if len(elem.transform) <= maxlevel+1:
            domain = domain.refined_by((elem.transform,))

    return domain

def portional_marking(domain, indicators, num, maxlevel=10):

    assert indicators.type == 'elementbased' or indicators.type == 'functionbased'
    assert isinstance(num, int), 'num should be an integer' 

    values = [abs(val) for val in indicators.indicators.values()]
    keys   = [key for key in indicators.indicators.keys()]
    to_refine = numpy.argsort(values)[-num:]

    if indicators.type == 'elementbased':
        assert len(domain) >= num, 'Amount of elements to be refined should be lower than the total amount of elements'
        marked = numpy.array(domain.elements)[to_refine]
    elif indicators.type == 'functionbased':
        assert len(indicators.basis) >= num, 'Amount of functions to be refined should be lower than the total amount of functions'
        mask = numpy.zeros(len(indicators.basis), dtype=bool)
        mask[to_refine] = True
        marked = domain.supp(indicators.basis, mask)

    for elem in marked:
        if len(elem.transform) <= maxlevel+1:
            domain = domain.refined_by((elem.transform,))

    return domain



    #     assert num < 1, 'Fraction should be lower than 1'
    #     row = sorted(list(indicators.indicators.values()))
    #     total = sum(indicators.indicators.values())
    #     Sum = 0
    #     for i, val in enumerate(row[::-1]):
    #         Sum += val
    #         if Sum >= total*num:
    #            threshold = row[-i] 
    #            break
    #            
    #     #threshold = max(list(indicators.indicators.values()))*(1-num)

    # to_refine = []
    # 
    # for key in list(indicators.indicators.keys()):
    #     if indicators.indicators[key] >= threshold:
    #         to_refine += [key]
    #         
    # if indicators.type == 'elementbased':
    #     assert len(domain) == len(indicators.indicators), 'Amount of indicators and amount of elements are not equal'
    #     marked = numpy.array(domain.elements)[to_refine]
    # elif indicators.type == 'functionbased':
    #     assert len(indicators.basis) == len(indicators.indicators), 'Amount of indicators and amount of basis functions are not equal'
    #     mask = numpy.zeros(len(indicators.basis), dtype=bool)
    #     mask[to_refine] = True
    #     marked = domain.supp(indicators.basis, mask)

    # for elem in marked:
    #     if len(elem.transform) <= maxlevel+1:
    #         domain = domain.refined_by((elem.transform,))

