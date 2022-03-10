# cython: language_level=3, boundscheck=False
"""
Copyright (c) 2015 Genome Research Ltd.

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

from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.map cimport map as mapcpp
from cython.operator cimport dereference as deref

cdef class WeightedChoice:
    def __cinit__(self):
        self.thisptr = new Chooser()
        self.pos = 0
    
    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr
    
    def __len__(self):
        return self.thisptr.len()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.pos >= len(self):
            self.pos = 0
            raise StopIteration
        
        iter = self.thisptr.iter(self.pos)
        self.pos += 1
        
        return {"pos": iter.pos, "ref": iter.ref.decode('utf8'),
            "alt": iter.alt.decode('utf8'), "prob": iter.prob, "offset": iter.offset}
    
    def append(self, WeightedChoice other):
        ''' combines the sites from two WeightedChoice objects
    
        NOTE: This is inefficient, but I can't figure out an easier way around
        NOTE: cython.
    
        Args:
            other: WeightedChoice object
        '''
    
        self.thisptr.append(deref(other.thisptr))
	
    def set_inherited_choices(self, vector[int] choices):
        self.thisptr.set_inherited_choices(choices)
	
    def add_choice(self, site, prob, ref='N', alt='N', offset=0, score = 1):
        """ add another possible choice for selection
        
        Args:
            site: a CDS position of the selected base.
            prob: probability of selecting this base.
            ref: string for reference allele
            alt: string for alternate allele
        """
        
        ref = ref.encode('utf8')
        alt = alt.encode('utf8')
        
        if len(ref) > 1 or len(alt) > 1:
            raise TypeError("requires single base alleles: {}, {}".format(ref, alt))
        
        self.thisptr.add_choice(site, prob, ref, alt, offset,score)

    def get_prob(self, site):
        return self.thisptr.get_prob(site)
	
    def choice(self):
        """ chooses a random element using a set of probability weights
        
        Returns:
            the name of the randomly selected element (e.g. position)
        """
        
        return self.choice_with_alleles()["pos"]
    
    def choice_with_alleles(self):
        """ chooses a random element, but include alleles in output
        
        Returns:
            a dictionary of the randomly selected element, with "pos", "ref"
            and "alt" entries.
        """
        
        choice = self.thisptr.choice()
        
        return {"pos": choice.pos, "ref": choice.ref.decode('utf8'),
            "alt": choice.alt.decode('utf8'), "offset": choice.offset, "prob": choice.prob}
    
    def get_summed_rate(self):
        """ return the cumulative probability for the object
        """
        return self.thisptr.get_summed_rate()

    def get_inherited_variants(self):
        """ return the control variants
        """
        return self.thisptr.get_inherited_variants()

    def get_inherited_count(self):
        return self.thisptr.get_inherited_count()


cdef extern from "simulate.h":
    void _get_distances(vector[int], vector[double] &)
    void _get_distances(vector[int], vector[vector[double]], vector[double] &)
    void _get_distances_coevol(vector[int], vector[double], vector[double] &)    
    #void _get_distances(vector[int], double[:,:], vector[double] &)
    vector[double] _scale_distances(vector[double], vector[double])
#    bool _has_zero(vector[int])
    bool _has_zero(vector[double])
    double _geomean(vector[double], double)
    double _geomean(vector[double], double, double)        
    bool _halt_permutation(double, int, double, double)
    vector[double] _simulate_distribution(Chooser, int, int, double, string)
    double _analyse_de_novos(Chooser, int, int, double, double, string)
    double _analyse_de_novos_entropy(Chooser, int, int, double, double, string)
    double _analyse_de_novos(Chooser, vector[vector[double]], int, int, double, double, string)
    double _analyse_de_novos_multi(vector[string],vector[string], vector[int], mapcpp[string,Chooser], vector[mapcpp[int,vector[double]]], int, int, double, double, string)
    double _analyse_de_novos_west(Chooser, vector[vector[double]], int, double, int, double, double)    
    double _analyse_de_novos_coevol(Chooser, vector[double], int, int, double, double, string)    
    #double _analyse_de_novos(Chooser, double[:,:], int, int, double)

def get_distances_1d(vector[int] positions):
    """ gets the distances between two or more CDS positions
    
    Args:
        positions: list of CDS positions as ints
    
    Returns:
        list of pairwise distances between sites
    """
    cdef vector[double] distances
    _get_distances(positions, distances)
    return distances
    
def get_distances(vector[int] positions, list three_d_locations):

    #cdef mapcpp[int,vector[double]] _three_d_locations = {int(k) : v for k,v in three_d_locations.iteritems()}
    cdef vector[vector[double]] _three_d_locations = three_d_locations   
    #cdef double[:,:] _three_d_locations = three_d_locations
    cdef vector[double] distances
    _get_distances(positions, _three_d_locations, distances)
    return distances

def get_distances_coevol(vector[int] positions, list coevol):

    #cdef mapcpp[int,vector[double]] _coevol = {int(k) : v for k,v in coevol.iteritems()}
    cdef vector[double] _coevol = coevol   
    #cdef double[:,:] _coevol = coevol
    cdef vector[double] distances
    _get_distances_coevol(positions, _coevol, distances)
    return distances


def get_distances_old(vector[int] positions):
    """ gets the distances between two or more CDS positions
    
    Args:
        positions: list of CDS positions as ints
    
    Returns:
        list of pairwise distances between sites
    """
    cdef vector[double] distances
    _get_distances(positions, distances)
    return distances

def scale_distances(vector[double] distances, vector[double] scores):
    cdef vector[double] scaled_distances
    scaled_distances = _scale_distances(distances, scores)
    return scaled_distances

def has_zero(vector[double] distances):
    """ figure out whether any of the pairwise distances is zero
    """
    
    return _has_zero(distances)

def geomean(vector[double] distances, double p):
    """ gets the geometric mean distance between two or more CDS positions
    
    Args:
        distances: list of distances between CDS positions
    
    Returns:
        provides the mean distance of the pairwise distances
    """
    
    return _geomean(distances, p, 1)

def geomean(vector[double] distances, double p, double offset):
    """ gets the geometric mean distance between two or more CDS positions
    
    Args:
        distances: list of distances between CDS positions
    
    Returns:
        provides the mean distance of the pairwise distances
    """
    
    return _geomean(distances, p, offset)

def simulate_distribution(WeightedChoice choices, int iterations, int de_novos_count, double p, string dist_file_output):
    """
    """
    
    return _simulate_distribution(deref(choices.thisptr), iterations, de_novos_count, p, dist_file_output)

#def analyse_de_novos(WeightedChoice choices, dict three_d_locations, int iterations, int N, double observed_value):
def analyse_de_novos(WeightedChoice choices, list three_d_locations, int iterations, int N, double observed_value, double p,  string dist_file_output):
    #cdef mapcpp[int,vector[double]] _three_d_locations = {int(k) : v for k,v in three_d_locations.iteritems()}	
    cdef vector[vector[double]] _three_d_locations = three_d_locations
    #cdef double[:,:] _three_d_locations = three_d_locations
    print("degree value is" + str(p))
    return _analyse_de_novos(deref(choices.thisptr), _three_d_locations, iterations,  N, observed_value, p, dist_file_output)

def analyse_de_novos_west(WeightedChoice choices, list three_d_locations, int iterations, double lambda_, int N, double observed_value, double p):
    #cdef mapcpp[int,vector[double]] _three_d_locations = {int(k) : v for k,v in three_d_locations.iteritems()}	
    cdef vector[vector[double]] _three_d_locations = three_d_locations
    #cdef double[:,:] _three_d_locations = three_d_locations
    print("degree value is" + str(p))
    return _analyse_de_novos_west(deref(choices.thisptr), _three_d_locations, iterations,  lambda_, N, observed_value, p)

def analyse_de_novos_coevol(WeightedChoice choices, list coevol, int iterations, int N, double observed_value, double p, string dist_file_output):
    cdef vector[double] _coevol = coevol
    #cdef double[:,:] _coevol = coevol
    return _analyse_de_novos_coevol(deref(choices.thisptr), _coevol, iterations,  N, observed_value, p, dist_file_output)

def analyse_de_novos_1d(WeightedChoice choices, int iterations, int de_novos_count, double observed_value, double p,  string dist_file_output):
    """
    """
    
    return _analyse_de_novos(deref(choices.thisptr), iterations, de_novos_count, observed_value, p, dist_file_output)

def analyse_de_novos_entropy(WeightedChoice choices, int iterations, int de_novos_count, double observed_value, double p,  string dist_file_output):
    """
    """
    
    return _analyse_de_novos_entropy(deref(choices.thisptr), iterations, de_novos_count, observed_value, p, dist_file_output)


def analyse_de_novos_multi(list chains,
                           list proteins,
                           list de_novo_count_per_chain,
                           list choices,
                           list three_d_list_of_dict,
                           int iterations,
                           int de_novos_count,
                           double observed_value,
                           double p,
			   string dist_file_output):
    #chains = [a.decode('utf8') for a in chains]
    cdef vector[Chooser] choices_
    cdef mapcpp[string,Chooser] choices_map
    for choice in choices:
        choices_.push_back(deref((<WeightedChoice>choice["rate"]).thisptr))
        choices_map[choice["protein"]] = (deref((<WeightedChoice>choice["rate"]).thisptr))
    cdef vector[mapcpp[int,vector[double]]] three_d_ = three_d_list_of_dict
    return _analyse_de_novos_multi(chains,
                                   proteins,
                             de_novo_count_per_chain,
                             choices_map,
                             three_d_,
                             iterations,
                             de_novos_count,
                             observed_value,
                             p,
			     dist_file_output)
			     
                             