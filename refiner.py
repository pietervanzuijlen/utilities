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

def dorfler_marking(domain, indicators, num, basis, maxlevel=10, select_type=None):

    assert num > 0 and num <= 1, 'the fraction should be between 0 and 1'
    assert select_type in [None, 'same_level', 'highest_supp','supp_only'], 'Invalid selection type given'

    if len(indicators) == len(domain):
        indtype = 'elementbased'
    elif len(indicators) == len(basis):
        indtype = 'functionbased'
    else:
        assert False, 'No indicator type has been found. Check if the right domain or basis is given'

    if indtype == 'elementbased':

        transforms = []
        for trans in domain.transforms:
            transforms.append(trans)

        values = numpy.array(abs(indicators))

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

        if select_type == 'same_level':
            marked = () 
            for ielem in refindi:
                marked = marked + tuple(select_same_level(domain, basis, ielem))
        elif select_type == 'highest_supp':
            marked = () 
            for ielem in refindi:
                marked = marked + tuple(select_highest_supp(domain, basis, ielem, indicators))
        elif select_type == 'supp_only':
            marked = tuple(select_supp_only(domain, basis, refindi))
        elif select_type != None:
            print('Selection type not recognised. No selection type is applied.') 
            marked = tuple(transforms[i] for i in refindi if len(transforms[i]) <= maxlevel+1)
        else:
            marked = tuple(transforms[i] for i in refindi if len(transforms[i]) <= maxlevel+1)

    elif indtype == 'functionbased':

        values = numpy.array(abs(indicators))

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

        elemind = basis.get_support(numpy.array(refindi))
        marked = domain.transforms[elemind]

    initial_size = len(domain)

    domain = domain.refined_by(marked)

    if len(domain) == initial_size:
        log.user('Nothing is refined')

    return domain

def fractional_marking(domain, indicators, num, basis, maxlevel=10, select_type=None):

    assert num > 0 and num <= 1, 'num should be a float between 0 and 1' 
    assert select_type in [None, 'same_level', 'highest_supp'], 'Invalid selection type given'

    if len(indicators) == len(domain):
        indtype = 'elementbased'
    elif len(indicators) == len(basis):
        indtype = 'functionbased'
    else:
        assert False, 'No indicator type has been found. Check if the right domain or basis is given'

    values = abs(indicators)
    threshold = max(values)*(1-num)
    refindi = []

    for i, val in enumerate(values):
        if val >= threshold:
            refindi += [i]

    if indtype == 'elementbased':
        if select_type == 'same_level':
            marked = () 
            for ielem in refindi:
                marked = marked + tuple(select_same_level(domain, basis, ielem))
        elif select_type == 'highest_supp':
            marked = () 
            for ielem in refindi:
                marked = marked + tuple(select_highest_supp(domain, basis, ielem, indicators))
        elif select_type == 'supp_only':
            marked = tuple(select_supp_only(domain, basis, refindi))
        elif select_type != None:
            print('Selection type not recognised. No selection type is applied.') 
            marked = [domain.transforms[i] for i in refindi]
        else:
            marked = [domain.transforms[i] for i in refindi]
    elif indtype == 'functionbased':
        elemind = basis.get_support(numpy.array(refindi))
        marked = domain.transforms[elemind]

    initial_size = len(domain)

    for trans in marked:
        if len(trans) <= maxlevel+1:
            domain = domain.refined_by((trans,))

    if len(domain) == initial_size:
        log.user('Nothing is refined')

    return domain

def portional_marking(domain, indicators, num, basis, maxlevel=10, select_type=None):

    assert isinstance(num, int), 'num should be an integer' 
    assert select_type in [None, 'same_level', 'highest_supp'], 'Invalid selection type given'

    if len(indicators) == len(domain):
        indtype = 'elementbased'
    elif len(indicators) == len(basis):
        indtype = 'functionbased'
    else:
        assert False, 'No indicator type has been found. Check if the right domain or basis is given'

    refindi = numpy.argsort(abs(indicators))[-num:]

    if indtype == 'elementbased':
        assert len(domain) >= num, 'Amount of elements to be refined should be lower than the total amount of elements'
        if select_type == 'same_level':
            marked = () 
            for ielem in refindi:
                marked = marked + tuple(select_same_level(domain, basis, ielem))
        elif select_type == 'highest_supp':
            marked = () 
            for ielem in refindi:
                marked = marked + tuple(select_highest_supp(domain, basis, ielem, indicators))
        elif select_type == 'supp_only':
            marked = tuple(select_supp_only(domain, basis, refindi))
        elif select_type != None:
            print('Selection type not recognised. No selection type is applied.') 
            marked = [domain.transforms[i] for i in refindi]
        else:
            marked = [domain.transforms[i] for i in refindi]
    elif indtype == 'functionbased':
        assert len(basis) >= num, 'Amount of functions to be refined should be lower than the total amount of functions'
        elemind = basis.get_support(numpy.array(refindi))
        marked = domain.transforms[elemind]
         

    for trans in marked:
        if len(trans) <= maxlevel+1:
            domain = domain.refined_by((trans,))

    return domain
    

#Element selection functions

def select_same_level(domain, basis, ielem):
    funcs = basis.get_dofs(ielem)
    irefelems = basis.get_support(funcs)
    ilvl = len(domain.transforms[ielem])
    transforms = []
    for irefelem in irefelems:
        trans = domain.transforms[irefelem] 
        if len(trans) <= ilvl:
            transforms.append(trans)

    return transforms 

def select_highest_supp(domain, basis, ielem, indicators):
    funcs = basis.get_dofs(ielem)
    ind_sum = 0 
    iref = []
    ilvl = len(domain.transforms[ielem])
    for func in funcs:
        irefelems = basis.get_support(func)
        if sum(indicators[irefelems]) > ind_sum:
            ind_sum = sum(indicators[irefelems])
            iref = irefelems
    print(type(irefelems))
    unsorted = domain.transforms[irefelems]
    transforms = []
    for trans in unsorted:
        if len(trans) <= ilvl:
            transforms.append(trans)

    return transforms 

def select_supp_only(domain, basis, ielems):
    funcs = basis.get_dofs(ielems)
    irefelems = []
    for func in funcs:
        supp = basis.get_support(func)
        if all(i in ielems for i in supp):
            irefelems += [ielem for ielem in supp]
    transforms = [domain.transforms[i] for i in irefelems]

    return transforms 
