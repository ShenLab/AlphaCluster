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

import gzip

from scipy.stats import chi2
from intervaltree import IntervalTree

from denovonear.load_gene import construct_gene_object

def parse_header(header):
    header = header.strip().split('\t')
    
    indices = {}
    indices['transcript'] = header.index('transcript')
    indices['symbol'] = header.index('gene')
    indices['chrom'] = header.index('chr')
    indices['amino_acids'] = header.index('amino_acids')
    indices['ratio'] = header.index('obs_exp')
    indices['chisq'] = header.index('chisq_diff_null')
    
    return indices

def parse(line, idx):
    line = line.strip().split('\t')
    
    transcript = line[idx['transcript']].split('.')[0]
    symbol = line[idx['symbol']]
    chrom = line[idx['chrom']]
    pos = line[idx['amino_acids']]
    ratio = float(line[idx['ratio']])
    chisq = float(line[idx['chisq']])
    
    return transcript, symbol, chrom, pos, ratio, chisq

def load_regional_constraint(path):
    ''' load in constraint regions, and determine chromosome coordinates
    '''
    
    constraint = {}
    with gzip.open(path, 'rt') as handle:
        header = handle.readline()
        idx = parse_header(header)
        
        for line in handle:
            tx, symbol, chrom, pos, ratio, chisq = parse(line, idx)
            
            if symbol not in constraint:
                constraint[symbol] = {'tx': tx, 'chrom': chrom, 'regions': []}
            
            constraint[symbol]['regions'].append({'pos': pos, 'ratio': ratio,
                'chisq': chisq})
    
    return constraint

def aa_to_chrom(tx, region):
    ''' convert an amino acid region of a transcript to chromosomal coordinates
    
    Args:
        tx: Transcript object for a gene
        region: start and end amino acid positions (dash-separated) e.g. '1-260'
    
    Returns:
        tuple of start and end chromosomal coordinates
    '''
    start, end = region.split('-')
    start = (int(start) - 1) * 3
    end = ((int(end) - 1) * 3) + 2
    
    return tx.get_position_on_chrom(start), tx.get_position_on_chrom(end)

def get_constrained_positions(ensembl, constraint, symbol, threshold=1e-4, ratio_threshold=1.0):
    ''' get the positions in the constrained regions
    '''
    
    regions = IntervalTree()
    
    if symbol not in constraint:
        return regions
    
    data = constraint[symbol]
    tx = construct_gene_object(ensembl, data['tx'])
    
    for region in data['regions']:
        p_value = chi2.sf(region['chisq'], df=1)
        if p_value > threshold:
            continue
        
        if region['ratio'] > ratio_threshold:
            continue
        
        start, end = aa_to_chrom(tx, region['pos'])
        if tx.get_strand() == '-':
            start, end = end, start
        
        regions[start:end + 1] = True
    
    return regions
