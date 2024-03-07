#!/usr/bin/env python3

import sys
import os
import itertools
import numpy as np

from .libs.pdb import Atom, Nucleotide, FROM_ANGSTROM_TO_OXDNA
from .libs import base

def print_usage():
    print("USAGE:", file=sys.stderr)
    print("\t%s PDB_file direction [use 35 or 53 if the PDB file has the nucleotides listed in the 3' -> 5' or 5' -> 3' direction, respectively]" % sys.argv[0], file=sys.stderr)
    print("\t[-m/--models-as-strands]", file=sys.stderr)
    exit(1)
    
def parse_options():
    shortArgs = 'm'
    longArgs = ['models-as-strands',]
    
    opts = {
        "PDB_file" : "",
        "oxDNA_direction" : True,
        "models_as_strands" : False,
    }
    
    try:
        import getopt
        args, positional_args = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-m' or k[0] == '--models-as-strands':
                print("## Different models will be interpreted as different strands", file=sys.stderr)
                opts["models_as_strands"] = True
            
        opts['PDB_file'] = positional_args[0]
        direction = positional_args[1]
        
        if direction == "35":
            opts["oxDNA_direction"] = True
        elif direction == "53":
            opts["oxDNA_direction"] = False
        else:
            print("The 'direction' argument should be either 35 or 53", file=sys.stderr)
            exit(1)
    except Exception:
        print_usage()
        
    return opts


def main():
    if len(sys.argv) < 3:
        print_usage()
        
    opts = parse_options()
        
    pdb_file = opts['PDB_file']
    oxDNA_direction = opts['oxDNA_direction']
    models_as_strands = opts['models_as_strands']
        
    pdb_strands = []
    with open(pdb_file) as f:
        strand = []
        old_residue = ""
        old_chain = ""
        for line in f.readlines():
            line = line.strip()
            if line.startswith("ATOM"):
                na = Atom(line)
                if old_chain != "":
                    if na.chain_id != old_chain and len(strand) != 0:
                        print("WARNING: a TER statement separating different strands (%s and %s) is missing" % (na.chain_id, old_chain), file=sys.stderr)
                        pdb_strands.append(strand)
                        strand = []
                    elif na.chain_id == old_chain and len(strand) == 0:
                        print("WARNING: a TER statement separates strands having the same chain id (%s)" % na.chain_id, file=sys.stderr)
                        
                if na.alternate != "":
                    if na.alternate == "A" or na.alternate == "1":
                        print("Alternate location for atom '%s' of residue '%s' encountered, using the line marked with the '%s' character." % (na.name, na.residue, na.alternate), file=sys.stderr)
                if na.residue_idx != old_residue:
                    nn = Nucleotide(na.residue, na.residue_idx)
                    if oxDNA_direction:
                        strand.append(nn)
                    else:
                        strand.insert(0, nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)
                old_chain = na.chain_id
            elif line.startswith("MODEL"):
                if not models_as_strands:
                    N_model = line.split()[1]
                    print("MODEL line detected: using the first MODEL encountered (%s)" % (N_model), file=sys.stderr)
            elif line.startswith("ENDMDL"):
                if not models_as_strands:
                    # by default we treat ENDMDL as the end of the file 
                    break;
                else:
                    # otherwise we treat it as the end of the strand
                    if len(strand) > 0:
                        pdb_strands.append(strand)
                        strand = []
            elif line.startswith("TER"):
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
    
    print("Using a box of size %g in oxDNA units (twice as big as the PDB bounding box size)" % (L), file=sys.stderr)
    
    system = base.System(box)
    strand = base.Strand()
    
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
    
    print("## Wrote data to '%s' / '%s'" % (configuration_file, topology_file), file=sys.stderr)
    print("## DONE", file=sys.stderr)
    

if __name__ == '__main__':
    main()