""" unit test the SequenceMethods class
"""

from __future__ import division
import unittest

from sequence_methods import SequenceMethods
from interval import Interval

class TestSequenceMthodsPy(unittest.TestCase):
    """ unit test the SequenceMethods class
    """
    
    def setUp(self):
        """ construct an interval object to add sequence to
        """
        
        chrom = "1"
        start = 100
        end = 200
        name = "TEST"
        strand = "+"
        exons = [(100, 120), (180, 200)]
        cds = [(110, 120), (180, 190)]
        
        self.gene = Interval(name, start, end, strand, chrom, exons, cds)
    
    def test_get_position_in_cds(self):
        """ test get_position_in_cds() works correctly
        """
        
        self.assertEqual(self.gene.get_position_in_cds(110), 0)
        self.assertEqual(self.gene.get_position_in_cds(120), 10)
        self.assertEqual(self.gene.get_position_in_cds(190), 21)
        
        self.gene.strand = "-"
        self.assertEqual(self.gene.get_position_in_cds(190), 0)
        self.assertEqual(self.gene.get_position_in_cds(180), 10)
        self.assertEqual(self.gene.get_position_in_cds(110), 21)
        
        self.gene.strand = "1"
        self.assertRaises(ValueError, self.gene.get_position_in_cds, 110)
    
    def test_get_codon_number_for_cds_position(self):
        """ test that get_codon_number_for_cds_position() works correctly
        """
        
        self.assertEqual(self.gene.get_codon_number_for_cds_position(0), 0)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(1), 0)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(2), 0)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(3), 1)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(4), 1)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(30000), 10000)
    
    def test_get_position_within_codon(self):
        """ test that get_position_within_codon() works correctly
        """
        
        self.assertEqual(self.gene.get_position_within_codon(0), 0)
        self.assertEqual(self.gene.get_position_within_codon(1), 1)
        self.assertEqual(self.gene.get_position_within_codon(2), 2)
        self.assertEqual(self.gene.get_position_within_codon(3), 0)
    
    def test_add_genomic_sequence(self):
        """ test that add_genomic_sequence() works correctly
        """
        
        self.gene.start = 0
        self.gene.end = 10
        self.gene.exons = [(0, 4), (6, 10)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_start = 2
        self.gene.cds_end = 8
        
        gdna = "AAAGGCCTTT"
        self.gene.cds_sequence = "AGGCTT"
        
        self.gene.add_genomic_sequence(gdna)
        self.assertEqual(self.gene.cds_sequence, "AGGCTT")
        
        # now check for a gene on the reverse strand
        self.gene.strand = "-"
        self.gene.cds_sequence = "AAGCCT"
        self.gene.add_genomic_sequence(gdna)
        self.assertEqual(self.gene.cds_sequence, "AAGCCT")
        
        
        
        
    
      
        
        
    

if __name__ == '__main__':
    unittest.main()

