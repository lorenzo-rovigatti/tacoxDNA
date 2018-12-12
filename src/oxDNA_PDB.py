#!/usr/bin/env python

import sys
import os
import numpy as np
import copy
from math import sqrt, sin

from libs.pdb import Atom, Nucleotide, FROM_OXDNA_TO_ANGSTROM
import libs.base as base
import libs.utils as utils
from libs.readers import LorenzoReader

DD12_PDB_PATH = "dd12.pdb"

def print_usage():
        print >> sys.stderr, "USAGE:"
        print >> sys.stderr, "\t%s topology configuration direction" % sys.argv[0]
        print >> sys.stderr, "\t[-H\--hydrogens=True] [-u\--uniform-residue-names]"
        exit(1)

def parse_options():
    shortArgs = 'H:u'
    longArgs = ['hydrogens=','uniform-residue-names']
    
    opts = {
        "configuration" : "",
        "topology" : "",
        "oxDNA_direction" : True,
        "print_hydrogens" : True,
        "uniform_residue_names" : False
    }
    
    try:
        import getopt
        args, positional_args = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-H' or k[0] == '--hydrogens':
                k_arg = k[1].lower()
                if k_arg.lower() == "true":
                    opts["print_hydrogens"] = True
                elif k_arg == "false":
                    print >> sys.stderr, "## Hydrogen atoms will *not* be printed"
                    opts["print_hydrogens"] = False
                else:
                    print >> sys.stderr, "The argument of '%s' should be either 'true' or 'false' (got '%s' instead)" % (k[0], k[1])
                    exit(1)
            elif k[0] == '-u' or k[0] == '--uniform-residue-names':
                    opts["uniform_residue_names"] = True
            
        opts['topology'] = positional_args[0]
        opts['configuration'] = positional_args[1]
        direction = positional_args[2]
        
        if direction == "35":
            opts["oxDNA_direction"] = True
        elif direction == "53":
            opts["oxDNA_direction"] = False
        else:
            print >> sys.stderr, "The 'direction' argument should be either 35 or 53"
            exit(1)
    except Exception:
        print_usage()
        
    return opts

def align(full_base, ox_base):
        theta = utils.get_angle(full_base.a3, ox_base._a3)
        # if the two bases are already essentially aligned then we do nothing
        if sin(theta) > 1e-3:
            axis = np.cross(full_base.a3, ox_base._a3)
            axis /= sqrt(np.dot(axis, axis))
            R = utils.get_rotation_matrix(axis, theta)
            full_base.rotate(R)
    
        theta = utils.get_angle(full_base.a1, ox_base._a1)
        if sin(theta) > 1e-3:
            axis = np.cross(full_base.a1, ox_base._a1)
            axis /= sqrt(np.dot(axis, axis))
            R = utils.get_rotation_matrix(axis, theta)
            full_base.rotate(R)

if __name__ == '__main__':
    options = parse_options()
        
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
        lr = LorenzoReader(options['topology'], options['configuration'])
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
    
    out_name = options['configuration'] + ".pdb"
    out = open (out_name, "w")
    
    print >> out, "HEADER    t=1.12"
    print >> out, "MODEL     1"
    print >> out, "REMARK ## 0,0"

    current_base_identifier = 'A'    
    for strand in s._strands:
        strand_pdb = []
        nucleotides_in_strand = strand._nucleotides
        if not options['oxDNA_direction']:
                nucleotides_in_strand = reversed(nucleotides_in_strand)
        for n_idx, nucleotide in enumerate(nucleotides_in_strand, 1):
            nb = base.number_to_base[nucleotide._base]
            my_base = copy.deepcopy(bases[nb])
            my_base.chain_id = s._nucleotide_to_strand[nucleotide.index]
            residue_suffix = ""
            if options["uniform_residue_names"] == False:
                # 3' end
                if nucleotide == strand._nucleotides[0] and not strand.is_circular(): 
                    residue_suffix = "3"
                # 5' end
                elif nucleotide == strand._nucleotides[-1]: 
                    residue_suffix = "5"

            my_base.idx = (nucleotide.index % 12) + 1
            align(my_base, nucleotide)
            my_base.set_base((nucleotide.pos_base - com) * FROM_OXDNA_TO_ANGSTROM)

            if correct_for_large_boxes:
                my_base.correct_for_large_boxes(box_angstrom)
            
            serial_residue = n_idx % 9999
            base_identifier = current_base_identifier
            nucleotide_pdb = my_base.to_pdb(base_identifier, options['print_hydrogens'], serial_residue, residue_suffix)
            strand_pdb.append(nucleotide_pdb)
            
                
        print >> out, "\n".join(x for x in strand_pdb)
        print >> out, "TER"
        
        if current_base_identifier == 'Z':
            current_base_identifier = 'A'
        else:
            current_base_identifier = chr(ord(current_base_identifier) + 1)
        
    print >> out, "ENDMDL"
    out.close()
    
    print >> sys.stderr, "## Wrote data to '%s'" % out_name
    print >> sys.stderr, "## DONE"
