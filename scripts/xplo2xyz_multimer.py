#/usr/bin/python
from optparse import OptionParser
import re
import os.path
import sys
 
class PDBAtom(object):
    def __init__(self, string):
        #this is what we need to parse
        #ATOM      1  CA  ORN     1       4.935   1.171   7.983  1.00  0.00      sega
        #XPLOR pdb files do not fully agree with the PDB conventions 
        self.name = string[12:16].strip()
        self.x = float(string[30:38].strip())
        self.y = float(string[38:46].strip())
        self.z = float(string[46:54].strip())
        self.residue_number = int(string[22:27].strip())
        self.residue = string[16:20].strip()
        self.chain =  string[21:22].strip()
        self.warnings = []
        if len(string) < 78:
            self.element = self.name[0]
            self.warnings.append('Chemical element name guessed ' +\
                                'to be %s from atom name %s' % (self.element, self.name))
        else:
            self.element = string[76:78].strip()
 
usage = 'usage: %prog [options] <file.pdb> [<file.multimer>]\n\n' \
        + '\tConvert protein data bank PDB file created by NIH-XPLOR to MULTIMER file\n' \
        + '\tReferences: http://www.wwpdb.org/docs.html, ' \
        + 'http://en.wikipedia.org/wiki/XYZ_file_format\n\n' \
        + '\tto get help type: python %prog -h'
version = '%prog 0.1 - convert pdb file to multimer'
opt = OptionParser(usage=usage,version=version)
opt.add_option('-o','--overwrite',action='store_true',dest='overwrite',\
                default=False, help='overwrite output file, if it exists')
opt.add_option('-v','--verbose',action='store_true', dest='verbose',\
                default=False, help='print info about files being processed')
opt.add_option('-c','--ca-only',action='store_true', dest='caonly',\
                default=True, help='if true, only print CA atoms')

(options, args) = opt.parse_args()
 
narg = len(args)
if narg == 0:
    opt.error('must provide name of pdb file')
elif narg > 2:
    opt.error('too many no-option arguments should be either one or two (second - name of multimer file)')
else:
    infile = args[0]
    pdb_re = re.compile('^(.+).pdb$', re.IGNORECASE)
    m = pdb_re.search(infile)
    if m:
        basename = m.group(1)
        if narg == 2:
            if args[1].endswith('.multimer'):
                outfile = args[1]
            else:
                opt.error('output file (second argument) must have .multimer extension - case insensitive')
        else:
            outfile = basename + '.multimer'
    else:
        opt.error('input file (first argument) must have .pdb extension - case insensitive')
 
if os.path.exists(outfile) and options.overwrite == False:
    opt.error('file %s exists, use -o or --overwrite otion to overwrite the file' % outfile)
 
if os.path.isfile(infile):
    pdb_file = open(infile,'r')
else:
    opt.error('file %s does not exist' % infile)
 
if options.verbose:
    sys.stderr.write('converting %s --> %s\n' % (infile, outfile))
 
lineno = 0
models = {}
atoms = []
#read pdb file
#start_of_new_protein = True
model_number = 0
for line in pdb_file:
    last_atom_number = 0
    lineno += 1
#    if line.startswith('ENDMODEL'):
#        start_of_new_protein = true
#        continue
#    if line.startswith('CONNECT'):
#        start_of_new_protein = false
#        continue
#    if line.startswith('MODEL') and start_of_new_protein:
#        models[model_number] = []
    if line.startswith('ATOM'):
        try:    
            atoms.append(PDBAtom(line))
            #models[model_number].append(PDBAtom(line))
        except:
            sys.stderr.write('\nProblem parsing line %d in file %s\n' % (lineno,infile))
            sys.stderr.write(line)
            sys.stderr.write('Probably ATOM entry is formatted incorrectly?\n')
            sys.stderr.write('Please refer to - http://www.wwpdb.org/documentation/format32/sect9.html#ATOM\n\n')
            sys.exit(1)
pdb_file.close()
 
#save multimer file
multimer_file = open(outfile,'w')
#multimer_file.write('%d\n' % len(atoms))
#multimer_file.write('multimer file converted from %s\n' % infile)
lineno = 2
num_hidden_warnings = 0
for atom in atoms:
    lineno += 1
    if options.caonly == True:
        if not atom.name == "CA":
            continue
        #print(atom.chain)
        #print(atom.residue_number)
        multimer_file.write('%s\t%d\t%f\t%f\t%f\n' % (atom.chain, atom.residue_number, atom.x, atom.y, atom.z, atom.residue))
    else:
        multimer_file.write('%s\t%f\t%f\t%f\n' % (atom.name, atom.x, atom.y, atom.z))
    if atom.warnings:
        if options.verbose:
            sys.stderr.write('Possible issue on line %d in %s\n' % (lineno, outfile))
            sys.stderr.write('\n'.join(atom.warnings))
            sys.stderr.write('\n')
        else:
            num_hidden_warnings += 1
 
multimer_file.close()
if options.verbose == False and num_hidden_warnings > 0:
    sys.stderr.write('file %s saved\n' % outfile)
    sys.stderr.write('%d warnings were not shown, ' % num_hidden_warnings)
    sys.stderr.write('please rerun with option -v to see them\n')
