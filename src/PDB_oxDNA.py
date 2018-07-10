#!/usr/bin/env python

import sys
import os
import itertools
import numpy as np

from libs.pdb import Atom, Nucleotide, FROM_ANGSTROM_TO_OXDNA
from libs import base

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >> sys.stderr, "Usage is %s PDB_file direction [use 35 or 53 if the PDB file has the nucleotides listed in the 3' -> 5' or 5' -> 3' direction, respectively]" % sys.argv[0]
        sys.exit(1)
        
    pdb_file = sys.argv[1]
    if sys.argv[2] == "35":
        oxDNA_direction = True
    elif sys.argv[2] == "53":
        oxDNA_direction = False
    else:
        print >> sys.stderr, "The second argument should be either 35 or 53"
        exit(1)
        
    pdb_strands = []
    with open(pdb_file) as f:
        strand = []
        old_residue = ""
        for line in f.readlines():
            line = line.strip()
            if line.startswith("ATOM"):
                na = Atom(line)
                if na.residue_idx != old_residue:
                    nn = Nucleotide(na.residue, na.residue_idx)
                    if oxDNA_direction:
                        strand.append(nn)
                    else:
                        strand.insert(0, nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)
            elif line.split()[0].strip() == "TER":
                pdb_strands.append(strand)
                strand = []
            # if the file does not contain any TER line we need to manually add the current strand to the list of strands
            elif line == "END" and len(pdb_strands) == 0 and len(strand) > 0:
                pdb_strands.append(strand)
                strand = []
    
    # sometimes files just end (without any END or TER line)
    if len(strand) > 0:
        pdb_strands.append(strand)
                
    box_low = np.array([1e6, 1e6, 1e6], dtype=np.float64)
    box_high = np.array([-1e6, -1e6, -1e6], dtype=np.float64)
    for nucl in itertools.chain(*pdb_strands):
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
    print len(pdb_strands)
    
    for pdb_strand in pdb_strands:
        strand = base.Strand()
        
        for nucl in pdb_strand:
            nucl.compute_as()
            
            com = nucl.get_com() * FROM_ANGSTROM_TO_OXDNA
            new_oxDNA_nucl = base.Nucleotide(com, nucl.a1, nucl.a3, nucl.base[0])
            strand.add_nucleotide(new_oxDNA_nucl)
            
        system.add_strand(strand, check_overlap=False)
                
    basename = os.path.basename(pdb_file)
    topology_file = basename + ".top"
    configuration_file = basename + ".oxdna"
    system.print_lorenzo_output(configuration_file, topology_file)
    
    print >> sys.stderr, "## Wrote data to '%s' / '%s'" % (configuration_file, topology_file)
    print >> sys.stderr, "## DONE"
    