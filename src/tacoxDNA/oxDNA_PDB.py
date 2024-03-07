#!/usr/bin/env python3

import sys
import os
import numpy as np
import copy
from math import sqrt, sin
from collections import defaultdict

from .libs.pdb import Atom, Nucleotide, FROM_OXDNA_TO_ANGSTROM
from .libs import base
from .libs import utils
from .libs.readers import LorenzoReader

DD12_PDB_PATH = "dd12_na.pdb"

def print_usage():
    print("USAGE:", file=sys.stderr)
    print("\t%s topology configuration direction" % sys.argv[0], file=sys.stderr)
    print("\t[-H\--hydrogens=True] [-u\--uniform-residue-names] [-o\--one-file-per-strand] [--rmsf-file]", file=sys.stderr)
    exit(1)

def parse_options():
    shortArgs = 'H:uo'
    longArgs = ['hydrogens=', 'uniform-residue-names', 'one-file-per-strand', 'rmsf-file=']

    opts = {
        "configuration" : "",
        "topology" : "",
        "oxDNA_direction" : True,
        "print_hydrogens" : True,
        "uniform_residue_names" : False,
        "one_file_per_strand" : False,
        "rmsf_bfactor" : "",
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
                    print("## Hydrogen atoms will *not* be printed", file=sys.stderr)
                    opts["print_hydrogens"] = False
                else:
                    print("The argument of '%s' should be either 'true' or 'false' (got '%s' instead)" % (k[0], k[1]), file=sys.stderr)
                    exit(1)
            elif k[0] == '-u' or k[0] == '--uniform-residue-names':
                opts["uniform_residue_names"] = True
            elif k[0] == '-o' or k[0] == '--one-file-per-strand':
                opts["one_file_per_strand"] = True
            elif k[0] == '--rmsf-file':
                opts["rmsf_bfactor"] = k[1]

        opts['topology'] = positional_args[0]
        opts['configuration'] = positional_args[1]
        direction = positional_args[2]

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


def main():
    opts = parse_options()

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
        else:
            bases[n.base] = n

    for n in nucleotides:
        n.a1, n.a2, n.a3 = utils.get_orthonormalized_base(n.a1, n.a2, n.a3)

    try:
        lr = LorenzoReader(opts['topology'], opts['configuration'])
        s = lr.get_system()
    except Exception as e:
        print("Parser error: %s" % e, file=sys.stderr)
        exit(1)

    rmsf_file = opts["rmsf_bfactor"]
    if rmsf_file:
        with open(rmsf_file) as f:
            try:
                # https://github.com/sulcgroup/oxdna_analysis_tools
                substrings = f.read().split("[")[1].split("]")[0].split(",")
            except Exception as e:
                print("Parsing error in RMSF file. Invalid Format: %s" % e, file=sys.stderr)
                exit(1)
            try:
                rmsf_per_nucleotide = {i: float(s)
                                       for i, s in enumerate(substrings)}
            except Exception as e:
                print("Parsing error in RMSF file. Conversion to float failed : %s" % e, file=sys.stderr)
                exit(1)
    else:
        rmsf_per_nucleotide = defaultdict(lambda: 1.00)

    ox_nucleotides = []
    s.map_nucleotides_to_strands()
    com = np.array([0., 0., 0.])
    for my_strand in s._strands:
        com += my_strand.cm_pos
    com /= s.N_strands

    box_angstrom = s._box * FROM_OXDNA_TO_ANGSTROM
    correct_for_large_boxes = False
    if np.any(box_angstrom[box_angstrom > 999]):
        print("At least one of the box sizes is larger than 999: all the atoms which are outside of the box will be brought back through periodic boundary conditions", file=sys.stderr)
        correct_for_large_boxes = True

    if opts['one_file_per_strand']:
        out_name = opts['configuration'] + "_1.pdb"
    else:
        out_name = opts['configuration'] + ".pdb"

    out = open(out_name, "w")

    current_base_identifier = 'A'
    for s_id, strand in enumerate(s._strands):
        strand_pdb = []
        nucleotides_in_strand = strand._nucleotides
        if not opts['oxDNA_direction']:
            nucleotides_in_strand = reversed(nucleotides_in_strand)
        for n_idx, nucleotide in enumerate(nucleotides_in_strand, 1):
            nb = base.number_to_base[nucleotide._base]
            my_base = copy.deepcopy(bases[nb])
            my_base.chain_id = s._nucleotide_to_strand[nucleotide.index]
            residue_type = ""

            if not strand.is_circular():
                # 3' end
                if nucleotide == strand._nucleotides[0]:
                    residue_type = "3"
                # 5' end
                elif nucleotide == strand._nucleotides[-1]:
                    residue_type = "5"

            if opts["uniform_residue_names"] == True:
                residue_suffix = ""
            else:
                residue_suffix = residue_type

            align(my_base, nucleotide)
            my_base.set_base((nucleotide.pos_base - com) * FROM_OXDNA_TO_ANGSTROM)

            if correct_for_large_boxes:
                my_base.correct_for_large_boxes(box_angstrom)

            residue_serial = n_idx % 9999
            base_identifier = current_base_identifier
            nucleotide_pdb = my_base.to_pdb(base_identifier, opts['print_hydrogens'], residue_serial, residue_suffix,
                                            residue_type, bfactor=rmsf_per_nucleotide[nucleotide.index])
            strand_pdb.append(nucleotide_pdb)


        print("\n".join(x for x in strand_pdb), file=out)
        print("TER", file=out)

        if opts['one_file_per_strand']:
            out.close()
            print("## Wrote strand %d's data to '%s'" % (s_id + 1, out_name), file=sys.stderr)
            # open a new file if needed
            if strand != s._strands[-1]:
                out_name = opts['configuration'] + "_%d.pdb" % (s_id + 2, )
                out = open(out_name, "w")
        else:
            # we update the base identifier only if a single file is printed
            if current_base_identifier == 'Z':
                current_base_identifier = 'A'
            else:
                current_base_identifier = chr(ord(current_base_identifier) + 1)

    if not opts['one_file_per_strand']:
        out.close()
        print("## Wrote data to '%s'" % out_name, file=sys.stderr)

    print("## DONE", file=sys.stderr)


if __name__ == '__main__':
    main()
