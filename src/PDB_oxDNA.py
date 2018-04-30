#!/usr/bin/env python

# TODO: need to improve parsing to support the TER statement

import sys
import os
import numpy as np

from libs.pdb import Atom, Nucleotide, FROM_ANGSTROM_TO_OXDNA
from libs import base

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >> sys.stderr, "Usage is %s PDB_file direction [use 35 or 53 if the PDB file has the nucleotides listed in the 3' -> 5' or 5' -> 3' direction, respectively]" % sys.argv[0]
        sys.exit(1)
        
    pdb_file = sys.argv[1]
    if sys.argv[2] != "35":
        oxDNA_direction = True
    elif sys.argv[2] != "53":
        oxDNA_direction = False
    else:
        print >> sys.stderr, "The second argument should be either 35 or 53"
        exit(1)
        
    with open(pdb_file) as f:
        nucleotides = []
        old_residue = ""
        for line in f.readlines():
            if line.startswith("ATOM"):
                na = Atom(line)
                if na.residue_idx != old_residue:
                    nn = Nucleotide(na.residue, na.residue_idx)
                    nucleotides.append(nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)
                
    # PDB files might list nucleotide in the 5' -> 3' order. oxDNA always uses the 3' -> 5' convention, hence the need to revert the array
    if not oxDNA_direction:
        nucleotides.reverse()
        
    box_low = np.array([1e6, 1e6, 1e6], dtype=np.float64)
    box_high = np.array([-1e6, -1e6, -1e6], dtype=np.float64)
    for nucl in nucleotides:
        com = nucl.get_com()
        for i in range(3):
            if com[i] < box_low[i]:
                box_low[i] = com[i]
            elif com[i] > box_high[i]:
                box_high[i] = com[i]
                
    L = 2 * np.max(box_high - box_low) * FROM_ANGSTROM_TO_OXDNA
    box = np.array([L, L, L])
    
    print >> sys.stderr, "Using a box of size %g in oxDNA units (twice as big as the PDB bounding box size)" % (L)
    
    system = base.System(box)
    strand = base.Strand()
    
    for nucl in nucleotides:
        nucl.compute_as()
        
        com = nucl.get_com() * FROM_ANGSTROM_TO_OXDNA
        new_oxDNA_nucl = base.Nucleotide(com, nucl.a1, nucl.a3, nucl.base[0])
        strand.add_nucleotide(new_oxDNA_nucl)
        
        if len(nucl.base) > 2:
            print >> sys.stderr, "ERROR: invalid PDB base type '%s'" % nucl.base
            exit(1)
        elif nucl.base[-1] == "3":
            if strand.N != 1:
                print >> sys.stderr, "ERROR: malformed PDB file. Trying to add a 3' nucleotide (name: %s, index: %s) to a non-empty strand" % (nucl.name, nucl.idx)
                exit(1)
        elif nucl.base[-1] == "5" or nucl == nucleotides[-1]: # TODO the rhs of the or will be removed in the future
            system.add_strand(strand, check_overlap=False)
            strand = base.Strand()
                
    basename = os.path.basename(pdb_file)
    system.print_lorenzo_output(basename + ".oxdna", basename + ".top")
            