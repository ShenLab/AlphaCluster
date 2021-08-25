
import asyncio
from pathlib import Path
import unittest
import tempfile

from denovonear.rate_limiter import RateLimiter
from denovonear.load_gene import construct_gene_object
from denovonear.gencode import Gencode, _parse_gtfline, _open_gencode

async def call(func, *args, **kwargs):
    ''' call ensembl rest API function
    '''
    async with RateLimiter(15) as ensembl:
        return await func(ensembl, *args, **kwargs)

def _run(func, *args, **kwargs):
    return asyncio.get_event_loop().run_until_complete(call(func, *args, **kwargs))

def write_gtf(path, lines):
    with open(path, 'wt') as output:
        for line in lines:
            output.write(line)

class TestGencode(unittest.TestCase):
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent.parent /  "data"
        self.gtf_path = self.folder / 'example.grch38.gtf'
        self.fasta_path = self.folder / 'example.grch38.fa'
    
    def test_gencode_opens(self):
        ''' test we can open a gencode object
        '''
        gencode = Gencode(self.gtf_path, self.fasta_path)
        self.assertEqual(len(gencode), 1)
        gene = gencode['OR4F5']
        with self.assertRaises(KeyError):
            gencode['ZZZZZZZ']
        
    def test_gencode_matches_ensembl(self):
        ''' test thaqt gencode data matches ensembl
        '''
        gencode = Gencode(self.gtf_path, self.fasta_path)
        for gencode_tx in gencode['OR4F5'].transcripts:
            # get the transcript ID (but trim the version number)
            tx_id = gencode_tx.get_name().split('.')[0]
            ensembl_tx = _run(construct_gene_object, tx_id)
            self.assertEqual(gencode_tx.get_exons(), ensembl_tx.get_exons())
            self.assertEqual(gencode_tx.get_cds(), ensembl_tx.get_cds())
            self.assertEqual(gencode_tx.get_cds_sequence(), ensembl_tx.get_cds_sequence())
    
    def test_parse_gtf_gene_line(self):
        ''' test we can parse a GTF line for a gene feature
        '''
        line = 'chr1\tHAVANA\tgene\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4";' \
            'gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";' \
            'level 2; havana_gene "OTTHUMG00000001094.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'gene', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'',
            'transcript_type': b'',
            'is_principal': False,
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_transcript_line(self):
        ''' test we can parse a GTF line for a transcript feature
        '''
        line = 'chr1\tHAVANA\ttranscript\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; '\
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'transcript', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_principal': True,
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_exon_line(self):
        '''test we can parse a GTF line for an exon feature
        '''
        line = 'chr1\tHAVANA\texon\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; ' \
            'exon_id "ENSE00002319515.1"; level 2; protein_id "ENSP00000334393.3"; ' \
            'tag "basic"; transcript_support_level "NA"; tag "appris_principal_1"; ' \
            'tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.2"; ' \
            'havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'exon', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_principal': False,  ## exons don't get checked for principal tag
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_cds_line(self):
        '''test we can parse a GTF line for a CDS feature
        '''
        line = 'chr1\tHAVANA\tCDS\t69091\t70005\t.\t+\t0\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; ' \
            'exon_id "ENSE00002319515.1"; level 2; protein_id "ENSP00000334393.3"; ' \
            'tag "basic"; transcript_support_level "NA"; tag "appris_principal_1"; ' \
            'tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.2"; ' \
            'havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'CDS', 
            'start': 69091,
            'end': 70005,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_principal': False,  ## CDS don't get checked for principal tag
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_start_codon_line(self):
        '''test we can parse a GTF line for a start codon feature
        '''
        line = 'chr1\tHAVANA\tstart_codon\t69091\t69093\t.\t+\t0\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; ' \
            'gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; ' \
            'transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'start_codon', 
            'start': 69091,
            'end': 69093,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_principal': False,
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_stop_codon_line(self):
        '''test we can parse a GTF line for a stop codon feature
        '''
        line = 'chr1\tHAVANA\tstop_codon\t70006\t70008\t.\t+\t0\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN";' \
            'gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; ' \
            'transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'stop_codon', 
            'start': 70006,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_principal': False,
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_UTR_line(self):
        '''test we can parse a GTF line for a UTR feature
        '''
        line = 'chr1\tHAVANA\tUTR\t70006\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; ' \
            'gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; ' \
            'transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'UTR', 
            'start': 70006,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_principal': False,
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_minus_strand(self):
        '''test we can parse a GTF line for a feature on the minus strand
        '''
        line = 'chr1\tHAVANA\tCDS\t70006\t70008\t.\t-\t.\tgene_name "TEST";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'CDS', 
            'start': 70006,
            'end': 70008,
            'strand': b'-',
            'symbol': b'TEST',
            'tx_id': b'',
            'transcript_type': b'',
            'is_principal': False,
            }
        self.assertEqual(obj, expected)
    
    def test__open_gencode_multi_gene(self):
        '''test we can parse a GTF with multiple genes
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST1";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST1"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST1"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST1"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tgene\t10\t30\t.\t-\t.\tgene_name "TEST2";\n' \
                'chr1\tHAVANA\ttranscript\t10\t30\t.\t-\t.\ttranscript_id "ENST_B";gene_name "TEST2"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t30\t.\t-\t.\ttranscript_id "ENST_B" gene_name "TEST2"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tCDS\t15\t30\t.\t-\t.\ttranscript_id "ENST_B" gene_name "TEST2"; transcript_type "protein_coding"\n'
        
        with tempfile.NamedTemporaryFile() as temp:
            write_gtf(temp.name, lines)
            data = _open_gencode(temp.name)
        
        self.assertEqual(len(data), 2)
        symbol1, tx1, is_principal = data[0]
        symbol2, tx2, is_principal = data[1]
        self.assertNotEqual(symbol1, symbol2)
    
    def test__open_gencode_not_coding(self):
        '''test we can parse GTFs without protein coding transcripts
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST"; transcript_type "processed_transcript"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "processed_transcript"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "processed_transcript"\n' \
                        
        with tempfile.NamedTemporaryFile() as temp:
            write_gtf(temp.name, lines)
            data = _open_gencode(temp.name)
        
        self.assertEqual(len(data), 0)
        
        # but if we allow all transcript types (not just protein_coding, we)
        with tempfile.NamedTemporaryFile() as temp:
            write_gtf(temp.name, lines)
            data = _open_gencode(temp.name, coding_only=False)
        
        self.assertEqual(len(data), 1)
        
    def test__open_gencode_multi_transcript(self):
        '''test we can parse a GTF with multiple transcripts for the same gene
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tgene\t10\t30\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t30\t.\t-\t.\ttranscript_id "ENST_B";gene_name "TEST"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t30\t.\t-\t.\ttranscript_id "ENST_B" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tCDS\t15\t30\t.\t-\t.\ttranscript_id "ENST_B" gene_name "TEST"; transcript_type "protein_coding"\n'
        
        with tempfile.NamedTemporaryFile() as temp:
            write_gtf(temp.name, lines)
            data = _open_gencode(temp.name)
        
        self.assertEqual(len(data), 2)
        symbol1, tx1, is_principal = data[0]
        symbol2, tx2, is_principal = data[1]
        self.assertEqual(symbol1, symbol2)
        self.assertEqual(tx1.get_name(), 'ENST_A')
        self.assertEqual(tx1.get_cds(), [{'start': 15, 'end': 20}])
        self.assertEqual(tx2.get_name(), 'ENST_B')
        self.assertEqual(tx2.get_cds(), [{'start': 15, 'end': 30}])
        
    def test__open_gencode_multi_exon(self):
        '''test we can parse a GTF into transcripts
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t100\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t100\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\tUTR\t10\t15\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\texon\t30\t40\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tCDS\t30\t40\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\texon\t90\t100\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n' \
                'chr1\tHAVANA\tUTR\t90\t100\t.\t-\t.\ttranscript_id "ENST_A" gene_name "TEST"; transcript_type "protein_coding"\n'
        
        with tempfile.NamedTemporaryFile() as temp:
            write_gtf(temp.name, lines)
            data = _open_gencode(temp.name)
        
        self.assertEqual(len(data), 1)
        symbol, tx, is_principal = data[0]
        self.assertEqual(symbol, 'TEST')
        self.assertEqual(tx.get_name(), 'ENST_A')
        self.assertEqual(tx.get_strand(), '-')
        self.assertEqual(tx.get_exons(), [{'start': 10, 'end': 20}, {'start': 30, 'end': 40}, {'start': 90, 'end': 100}])
        self.assertEqual(tx.get_cds(), [{'start': 15, 'end': 20}, {'start': 30, 'end': 40}])
    

