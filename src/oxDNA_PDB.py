#!/usr/bin/env python

import sys
import os
import numpy as np
import copy
from math import sqrt

from libs.pdb import Atom, Nucleotide, FROM_OXDNA_TO_ANGSTROM
import libs.base as base
import libs.utils as utils
from libs.readers import LorenzoReader

DD12_PDB_PATH = "dd12.pdb"

def align(full_base, ox_base):
        theta = utils.get_angle(full_base.a3, ox_base._a3)
        axis = np.cross(full_base.a3, ox_base._a3)
        axis /= sqrt(np.dot(axis, axis))
        R = utils.get_rotation_matrix(axis, theta)
        full_base.rotate(R)
    
        theta = utils.get_angle(full_base.a1, ox_base._a1)
        axis = np.cross(full_base.a1, ox_base._a1)
        axis /= sqrt(np.dot(axis, axis))
        R = utils.get_rotation_matrix(axis, theta)
        full_base.rotate(R)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >> sys.stderr, "Usage is %s topology configuration" % sys.argv[0]
        sys.exit(1)

    with open(os.path.join(os.path.dirname(__file__), DD12_PDB_PATH)) as f:
        nucleotides = []
        old_residue = ""
        for line in f.readlines():
            if len(line) > 77:
                na = Atom(line)
                if na.residue_idx != old_residue:
                    nn = Nucleotide(na.residue, na.residue_idx)
                    nucleotides.append(nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)
           
    bases = {}
    for n in nucleotides:
        n.compute_as()
        if n.base in bases:
            if n.check < bases[n.base].check: bases[n.base] = copy.deepcopy(n)
        else: bases[n.base] = n
    
    for n in nucleotides:
        n.a1, n.a2, n.a3 = utils.get_orthonormalized_base(n.a1, n.a2, n.a3)
    
    try:
        lr = LorenzoReader(sys.argv[1], sys.argv[2])
        s = lr.get_system()
    except Exception as e:
        print >> sys.stderr, "Parser error: %s" % e
        exit(1)
    
    ox_nucleotides = []
    s.map_nucleotides_to_strands()
    com = np.array([0., 0., 0.])
    for my_strand in s._strands:
        com += my_strand.cm_pos
    com /= s.N_strands
    
    box_angstrom = s._box * FROM_OXDNA_TO_ANGSTROM
    correct_for_large_boxes = False
    if np.any(box_angstrom[box_angstrom > 999]):
        print >> sys.stderr, "At least one of the box sizes is larger than 999: all the atoms which are outside of the box will be brought back"
        correct_for_large_boxes = True
    
    out_name = sys.argv[2] + ".pdb"
    out = open (out_name, "w")
    
    print >> out, "HEADER    t=1.12"
    print >> out, "MODEL     1"
    print >> out, "REMARK ## 0,0"
    first_id = True
    old_chain_id = -1
    for nucleotide in s._nucleotides:
        nb = base.number_to_base[nucleotide._base]
        my_base = copy.deepcopy(bases[nb])
        my_base.chain_id = s._nucleotide_to_strand[nucleotide.index]
        is_3_prime = False
        if my_base.chain_id != old_chain_id:
            if old_chain_id != -1:
                print >> out, "TER"
            old_chain_id = my_base.chain_id
            is_3_prime = True
        my_base.idx = (nucleotide.index % 12) + 1
        align(my_base, nucleotide)
        my_base.set_base((nucleotide.pos_base - com) * FROM_OXDNA_TO_ANGSTROM)
        ox_nucleotides.append(my_base)
        if correct_for_large_boxes:
            my_base.correct_for_large_boxes(box_angstrom)
        print >> out, my_base.to_pdb("A", False, is_3_prime)
    print >> out, "REMARK ## 0,0"
    print >> out, "TER"
    print >> out, "ENDMDL"
    
    out.close()
    
    print >> sys.stderr, "## Wrote data to '%s'" % out_name
    print >> sys.stderr, "## DONE"
