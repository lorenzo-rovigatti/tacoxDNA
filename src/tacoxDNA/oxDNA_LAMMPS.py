#!/usr/bin/env python3

import numpy as np
import sys
from .libs.readers import LorenzoReader
from .libs.constants import mass_in_lammps, inertia_in_lammps, number_oxdna_to_lammps

def exyz_to_quat(mya1, mya3):
    mya2 = np.cross(mya3, mya1)
    myquat = [1, 0, 0, 0]

    q0sq = 0.25 * (mya1[0] + mya2[1] + mya3[2] + 1.0)
    q1sq = q0sq - 0.5 * (mya2[1] + mya3[2])
    q2sq = q0sq - 0.5 * (mya1[0] + mya3[2])
    q3sq = q0sq - 0.5 * (mya1[0] + mya2[1])

    # some component must be greater than 1/4 since they sum to 1
    # compute other components from it

    if q0sq >= 0.25:
        myquat[0] = np.sqrt(q0sq)
        myquat[1] = (mya2[2] - mya3[1]) / (4.0 * myquat[0])
        myquat[2] = (mya3[0] - mya1[2]) / (4.0 * myquat[0])
        myquat[3] = (mya1[1] - mya2[0]) / (4.0 * myquat[0])
    elif q1sq >= 0.25:
        myquat[1] = np.sqrt(q1sq)
        myquat[0] = (mya2[2] - mya3[1]) / (4.0 * myquat[1])
        myquat[2] = (mya2[0] + mya1[1]) / (4.0 * myquat[1])
        myquat[3] = (mya1[2] + mya3[0]) / (4.0 * myquat[1])
    elif q2sq >= 0.25:
        myquat[2] = np.sqrt(q2sq)
        myquat[0] = (mya3[0] - mya1[2]) / (4.0 * myquat[2])
        myquat[1] = (mya2[0] + mya1[1]) / (4.0 * myquat[2])
        myquat[3] = (mya3[1] + mya2[2]) / (4.0 * myquat[2])
    elif q3sq >= 0.25:
        myquat[3] = np.sqrt(q3sq)
        myquat[0] = (mya1[1] - mya2[0]) / (4.0 * myquat[3])
        myquat[1] = (mya3[0] + mya1[2]) / (4.0 * myquat[3])
        myquat[2] = (mya3[1] + mya2[2]) / (4.0 * myquat[3])

    norm = 1.0 / np.sqrt(myquat[0] * myquat[0] + myquat[1] * myquat[1] + \
			  myquat[2] * myquat[2] + myquat[3] * myquat[3])
    myquat[0] *= norm
    myquat[1] *= norm
    myquat[2] *= norm
    myquat[3] *= norm

    return np.array([myquat[0], myquat[1], myquat[2], myquat[3]])

def main():
    if len(sys.argv) < 3:
        print("Usage is %s topology configuration" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    try:
        lr = LorenzoReader(sys.argv[1], sys.argv[2])
        s = lr.get_system()
    except Exception as e:
        print("Parser error: %s" % e, file=sys.stderr)
        exit(1)
    
    s.map_nucleotides_to_strands()
    box = s._box
   
    # get total number of bonds
    N_bonds = 0
    for strand in s._strands:
        N_bonds += strand.get_lammps_N_of_bonds_strand()

    out_name = sys.argv[2] + ".lammps"
    out = open (out_name, "w")

    out.write('# LAMMPS data file\n')
    out.write('%d atoms\n' % s.N)
    out.write('%d ellipsoids\n' % s.N)
    out.write('%d bonds\n' % N_bonds)
    out.write('\n')
    out.write('4 atom types\n')
    out.write('1 bond types\n')
    out.write('\n')
    out.write('# System size\n')
    out.write('%f %f xlo xhi\n' % (0, box[0]))
    out.write('%f %f ylo yhi\n' % (0, box[1]))
    out.write('%f %f zlo zhi\n' % (0, box[2]))

    out.write('\n')
    out.write('Masses\n')
    out.write('\n')
    out.write('1 3.1575\n')
    out.write('2 3.1575\n')
    out.write('3 3.1575\n')
    out.write('4 3.1575\n')

    out.write('\n')
    out.write('# Atom-ID, type, position, molecule-ID, ellipsoid flag, density\n')
    out.write('Atoms\n')
    out.write('\n')

    for nucleotide in s._nucleotides:
        out.write('%d %d %22.15le %22.15le %22.15le %d 1 1\n' \
              % (nucleotide.index + 1, number_oxdna_to_lammps[nucleotide._base] + 1, \
                 nucleotide.cm_pos[0], nucleotide.cm_pos[1], nucleotide.cm_pos[2], \
                 s._nucleotide_to_strand[nucleotide.index] + 1))

    out.write('\n')
    out.write('# Atom-ID, translational, rotational velocity\n')
    out.write('Velocities\n')
    out.write('\n')

    for nucleotide in s._nucleotides:
        v_rescaled = np.array(nucleotide._v) / np.sqrt(mass_in_lammps)
        L_rescaled = np.array(nucleotide._L) * np.sqrt(inertia_in_lammps)
        out.write("%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n" \
              % (nucleotide.index + 1, v_rescaled[0], v_rescaled[1], v_rescaled[2], L_rescaled[0], L_rescaled[1], L_rescaled[2]))

    out.write('\n')
    out.write('# Atom-ID, shape, quaternion\n')
    out.write('Ellipsoids\n')
    out.write('\n')

    for nucleotide in s._nucleotides:
        quaternions = exyz_to_quat(nucleotide._a1, nucleotide._a3)
        out.write(\
        "%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n"  \
          % (nucleotide.index + 1, 1.1739845031423408, 1.1739845031423408, 1.1739845031423408, \
        quaternions[0], quaternions[1], quaternions[2], quaternions[3]))

    out.write('\n')
    out.write('# Bond topology\n')
    out.write('Bonds\n')
    out.write('\n')
    idx = 1
    for strand in s._strands:
        bonds = strand.get_lammps_bonds()
        for b in bonds:
            out.write("%d  %d  %s\n" % (idx, 1, b))
            idx += 1

    out.close()

    print("## Wrote data to '%s'" % out_name, file=sys.stderr)
    print("## DONE", file=sys.stderr)


if __name__ == '__main__':
    main()
    
