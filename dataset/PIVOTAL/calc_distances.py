import pandas as pd
import argparse
import sys
import os
from Bio.PDB import *
from Bio.PDB.vectors import *
import pickle


def calc_dihedral_atom(v1, v2, v3, v4):
    v1, v2, v3, v4 = Vector(v1.coord), Vector(v2.coord), Vector(
        v3.coord), Vector(v4.coord)
    return calc_dihedral(v1, v2, v3, v4)


def calc_angle_atom(v1, v2, v3):
    v1, v2, v3 = Vector(v1.coord), Vector(v2.coord), Vector(v3.coord)
    return calc_angle(v1, v2, v3)


def calc_geometry(r1, r2):
    if r1 == r2:
        return [0.0] * 4

    cb1 = r1['CA'] if 'CB' not in r1 else r1['CB']
    cb2 = r2['CA'] if 'CB' not in r2 else r2['CB']

    dis = cb1 - cb2

    try:
        omega = calc_dihedral_atom(r1['CA'], cb1, cb2, r2['CA'])
    except:
        omega = 0.0

    try:
        theta = calc_dihedral_atom(r1['N'], r1['CA'], cb1, cb2)
    except:
        theta = 0.0

    try:
        phi = calc_angle_atom(r1['CA'], cb1, cb2)
    except:
        phi = 0.0

    return [dis, omega, theta, phi]


def _get_res_id(r):
    het, resseq, icode = r.get_id()
    return (str(resseq) + icode).strip()


def calc_pdb(path, out):
    parser = PDBParser()
    struc = parser.get_structure('pdb', path)
    res = {}

    for a in struc.get_residues():
        aid = _get_res_id(a)
        for b in struc.get_residues():
            bid = _get_res_id(b)
            #if a.resname == 'GLY' or b.resname == 'GLY':
            #    continue
            res[(aid, bid)] = calc_geometry(a, b)

    with open(out, 'wb') as fw:
        pickle.dump(res, fw)


calc_pdb(sys.argv[1], sys.argv[2])
