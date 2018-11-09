from nutils import *


def amount(domain, indicators, N):

    assert len(domain) == len(indicators.indicators), 'Amount of indicators and amount of elements are not equal'
    assert len(domain) >= N, 'Amount of functions to be refined should be lower than the total amount of elements'

    to_refine = []
    max_values = sorted(list(indicators.indicators.values()))[-N:]

    for key in list(indicators.indicators.keys()):
        if indicators.indicators[key] in max_values:
            to_refine += [key]

    domain = domain.refined_by( e.transform for e in numpy.array(domain.elements)[to_refine] )

    return domain


def fraction(domain, indicators, frac):

    assert len(domain) == len(indicators.indicators), 'Amount of indicators and amount of elements are not equal'
    assert frac < 1, 'Fraction should be lower than 1'
    
    to_refine = []
    max_value = max(list(indicators.indicators.values()))

    for key in list(indicators.indicators.keys()):
        if indicators.indicators[key] >= max_value*(1-frac):
            to_refine += [key]

    domain = domain.refined_by( e.transform for e in numpy.array(domain.elements)[to_refine] )

    return domain

