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

def refine(domain, indicators, num, basis, maxlevel=10, grid=None, marker_type=None, select_type=None):

    indicators = abs(indicators)

    # Get indicater indices selected by the corresponding marker type
    if isinstance(num, int):
    # portional
        assert marker_type == None or marker_type =='portional', 'Invalid marker type or num has been giving'
        refindices = numpy.argsort(indicators)[-num:]
    elif isinstance(num, float) and (marker_type == None or marker_type == 'fractional'): 
    # fractional
        assert num > 0 and num < 1, 'Given num should be inbetween 0 and 1'
        threshold = max(indicators)*(1-num)
        refindices = numpy.array([i for i, val in enumerate(indicators) if val >= threshold])
    elif isinstance(num, float) and (marker_type == 'dorfler'): 
    # dorfler
        assert num > 0 and num < 1, 'Given num should be inbetween 0 and 1'
        threshold = indicators.sum()*num
        sortindi = numpy.flip(numpy.argsort(indicators))
        sortcumsum = indicators[sortindi].cumsum()
        i = numpy.argmax(sortcumsum>=threshold)
        refindices = sortindi[:i+1]
    else:
        assert False, 'Invalid marker type or num has been giving'

    # Get the transforms of elements to be refined
    if len(indicators) == len(domain):
    # elementbased
        if select_type==None:
            marked = domain.transforms[numpy.sort(refindices)] 
        elif select_type=='highest_supp':
            marked = select_highest_supp(domain, basis, refindices, indicators)
        elif select_type=='same_level':
            marked = select_same_level(domain, basis, refindices)
        else:
            assert False, 'Invalid selection type or num has been giving'

    elif len(indicators) == len(basis):
    # functionbased
        refelems = basis.get_support(refindices)
        marked = domain.transforms[refelems]
    else:
        assert False, 'No indicator type has been found. Check if the right domain or basis is given'

    # Discard elems which have reached maxlevel
    to_refine = tuple(trans for trans in marked if len(trans) <= maxlevel+1)
    if not to_refine:
        log.user('Nothing is refined')
        something_refined = False
    if to_refine:
        something_refined = True
        
    domain = domain.refined_by(to_refine)

    # Refinement with grid
    if grid:
        gridrefindices = numpy.array([grid.transforms.index(trans) for trans in to_refine])
        to_refine_grid = tuple(trans for trans in grid.transforms[gridrefindices])
        grid = grid.refined_by(tuple(to_refine_grid))

        return domain, grid, something_refined
    else: 
        return domain, something_refined


def select_same_level(domain, basis, ielems):
    marked = []
    for ielem in ielems:
        ilvl = len(domain.transforms[ielem])
        funcs = basis.get_dofs(ielem)
        irefelems = basis.get_support(funcs)
        marked += [trans for trans in domain.transforms[irefelems] if len(trans) <= ilvl]
    return marked 


def select_highest_supp(domain, basis, ielems, indicators):
    marked = []
    for ielem in ielems:
        ilvl = len(domain.transforms[ielem])
        funcs = basis.get_dofs(ielem)
        indsum = []
        for func in funcs:
            suppelems = basis.get_support(func)
            indsum += [sum(indicators[suppelems])]
        irefelems = basis.get_support(funcs[numpy.argmax(indsum)])
        marked += [trans for trans in domain.transforms[irefelems] if len(trans) <= ilvl]
    return marked 
