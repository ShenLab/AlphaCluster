"""
Copyright (c) 2017 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from itertools import groupby, count

from severity.open_mutations import LOF_CQ

def as_range(g):
    l = list(g)
    return l[0], l[-1]

def weight_site(site, consequence, scores, weights, constrained):
    ''' convert CADD score to the reweighted score
    
    Args:
        site: dict of data for a site. Includes 'pos' and 'alt' keys
        consequence: vep-style consequence for the site
        scores: dictionary of CADD scores, indexed by (position, alt) tuples
        weights: dict of weights for different consequence types and CADD
            thresholds.
        constrained: IntervalTree defining ranges within regional constraint for
            a gene.
    
    Returns:
        reweighted score
    '''
    
    score = scores[(site['pos'], site['alt'])]
    
    constraint = 'unconstrained'
    if site['pos'] in constrained:
        constraint = 'constrained'
    
    if consequence in LOF_CQ:
        return weights['truncating']
    else:
        for x in weights['altering'][constraint][score]:
            return x.data

def get_severity(cadd, chrom, rates, weights, constrained):
    ''' get CADD scores for a specific alt at a specific site
    
    See downloadable CADD files here: http://cadd.gs.washington.edu/download
    
    Args:
        cadd: pysam.TabixFile for quick fetching of CADD scores
        chrom: chromosome
        rates: Weighted Choice object for all sites in a gene
        weights:
        constrained: IntervalTree object for ranges under regional constraint
    
    Returns:
        list of weighted score at the given site for the given alt allele
    '''
    
    if type(rates) == dict:
        positions = sorted(set([ x['pos'] for cq in rates for x in rates[cq] ]))
    else:
        positions = sorted(set([ x['pos'] for x in rates ]))
    
    scores = {}
    for _, g in groupby(positions, key=lambda n, c=count(): n-next(c)):
        # define the start and end of a continuous set of positions, so we can
        # load cadd scores by iterating through the cadd file efficiently
        start, end = as_range(g)
        for line in cadd.fetch(chrom, start-1, end):
            _, pos, _, alt, _, score = line.split('\t')
            scores[(int(pos), alt)] = float(score)
    
    # match the cadd scores to the order of sites in the rates object
    if type(rates) == dict:
        cqs = {'synonymous': 'synonymous_variant', 'nonsense': 'stop_gained',
            'missense': 'missense_variant', 'splice_lof': 'splice_donor_variant'}
        return [ scores[x['pos'], x['alt']] for cq in sorted(rates) for x in rates[cq] ]
    else:
        return [ scores[x['pos'], x['alt']] for x in rates ]
    # # match the cadd scores to the order of sites in the rates object
    # if type(rates) == dict:
    #     cqs = {'synonymous': 'synonymous_variant', 'nonsense': 'stop_gained',
    #         'missense': 'missense_variant', 'splice_lof': 'splice_donor_variant'}
    #     return [ weight_site(x, cqs[cq], scores, weights, constrained) for cq in sorted(rates) for x in rates[cq] ]
    # else:
    #     return [ weight_site(x, x['consequence'], scores, weights, constrained) for x in rates ]