#!/usr/bin/env python3

import numpy as np
import sys, os
from .libs import base
from .libs import reader_lammps_init 
from .libs.constants import mass_in_lammps, inertia_in_lammps, number_oxdna_to_lammps


def quat_to_exyz(myquat):
    sqw = myquat[0] * myquat[0];
    sqx = myquat[1] * myquat[1];
    sqy = myquat[2] * myquat[2];
    sqz = myquat[3] * myquat[3];

    invs = 1 / (sqx + sqy + sqz + sqw)
    m00 = (sqx - sqy - sqz + sqw) * invs ;
    m22 = (-sqx - sqy + sqz + sqw) * invs ;
    
    tmp1 = myquat[1] * myquat[2];
    tmp2 = myquat[3] * myquat[0];
    m10 = 2.0 * (tmp1 + tmp2) * invs ;

    tmp1 = myquat[1] * myquat[3];
    tmp2 = myquat[2] * myquat[0];
    m20 = 2.0 * (tmp1 - tmp2) * invs ;
    m02 = 2.0 * (tmp1 + tmp2) * invs ;
    tmp1 = myquat[2] * myquat[3];
    tmp2 = myquat[1] * myquat[0];
    m12 = 2.0 * (tmp1 - tmp2) * invs ; 

    mya1 = np.array([m00, m10, m20])
    mya3 = np.array([m02, m12, m22])

    return mya1, mya3

def main():
    if len(sys.argv) < 2:
        print("USAGE:", file=sys.stderr)
        print("\t%s lammps_data_file [lammps_trajectory_file]" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    try:
        conf = reader_lammps_init.Lammps_parser(sys.argv[1])
    except ValueError as e:
        print(e, file=sys.stderr)
        exit(1)
        
    N = conf.natoms 
    box = np.array([0, 0., 0.])
    box[0] = conf.Lx
    box[1] = conf.Ly
    box[2] = conf.Lz

    system = base.System(box)

    strands = []
    for i in range(conf.nstrands):
        strands.append(base.Strand())

    for i in range(N):
        cm = conf.xyz[i,:]
        quaternions = conf.ellipsoids[i,:]
        a1, a3 = quat_to_exyz(quaternions)
        b = number_oxdna_to_lammps[(conf.bases[i]+3)%4] 

        v = np.array(conf.v[i,:]) * np.sqrt(mass_in_lammps)
        Lv = np.array(conf.Lv[i,:]) / np.sqrt(inertia_in_lammps)

        strands[conf.strand[i] - 1].add_nucleotide(base.Nucleotide(cm, a1, a3, b, b, v, Lv))

        # close strand 
        next_bond = conf.bonds[i][1]
        if next_bond != -1 and next_bond != i + 1:
            if conf.strand[i] != conf.strand[next_bond]:
                print("Wrong bond arising between two different strands", file=sys.stderr)
            else:
                strands[conf.strand[i] - 1].make_circular()

    for i in range(conf.nstrands):
        system.add_strand(strands[i])

    basename = os.path.basename(sys.argv[1])
    topology_file = basename + ".top"
    configuration_file = basename + ".oxdna"
    system.print_lorenzo_output(configuration_file, topology_file)

    # optional conversion of LAMMPS trajectory into native oxDNA format
    if len(sys.argv) == 3:
        with open(sys.argv[2], 'r') as lmptrj, open(configuration_file, 'w') as oxconf:
            line = lmptrj.readline()
            
            while line:
                if line.startswith('ITEM: TIMESTEP'):
                    t = int(lmptrj.readline())
            
                if line.startswith('ITEM: NUMBER OF ATOMS'):
                    natoms = int(lmptrj.readline())
                    if natoms != conf.natoms:
                        print("ERROR: A configuration stored in the trajectory file contains a number of nucleotides %d that differs from the %d in the datafile" % \
                          (natoms, conf.natoms), file=sys.stderr)            
                        sys.exit(1)
                if line.startswith('ITEM: BOX BOUNDS'):
                    line = lmptrj.readline()
                    xlo, xhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                    Lx = xhi - xlo
                    line = lmptrj.readline()
                    ylo, yhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                    Ly = yhi - ylo
                    line = lmptrj.readline()
                    zlo, zhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                    Lz = zhi - zlo
            
                if line.startswith('ITEM: ATOMS'):
                    aux = line.split()
            
                    # find column number in trajectory file
            
                    keyx = aux.index('x') - 2  # x
                    keyz = aux.index('z') - 1  # z (exclusive)
            
                    keyvx = aux.index('vx') - 2  # vx
                    keyvz = aux.index('vz') - 1  # vz (exclusive)
            
                    keylx = aux.index('angmomx') - 2  # angular momentum x
                    keylz = aux.index('angmomz') - 1  # angular momentum z (exclusive)
            
                    keyq0 = aux.index('c_quat[1]') - 2  # quat0
                    keyq3 = aux.index('c_quat[4]') - 1  # quat3 (exclusive)
            
                    N = natoms
            
                    # read position, velocity, quaternions, angular momentum
            
                    xyz = np.zeros((N, 3), dtype=float)
                    vel = np.zeros((N, 3), dtype=float)
                    quat = np.zeros((N, 4), dtype=float)
                    angmom = np.zeros((N, 3), dtype=float)
            
                    for n in range(N):
                        line = lmptrj.readline()
                        index = int(line.split()[0]) - 1
                        
                        xyz[index,:] = np.float32(line.split()[keyx:keyz])
                        vel[index,:] = np.float32(line.split()[keyvx:keyvz])
                        quat[index,:] = np.float32(line.split()[keyq0:keyq3])
                        angmom[index,:] = np.float32(line.split()[keylx:keylz]) 
            
                    # write oxDNA data to file
            
                    # header
                    oxconf.write('t = %d\n' % t)
                    oxconf.write('b = %f %f %f\n' % (Lx, Ly, Lz))
                    oxconf.write('E = 0.000000 0.000000 0.000000\n')
            
                    # atom data
                    for n in range(N):
                        cm = xyz[n,:]
                        quaternions = quat[n,:]
                        a1, a3 = quat_to_exyz(quaternions)
                        v = np.array(vel[n,:]) * np.sqrt(mass_in_lammps)
                        Lv = np.array(angmom[n,:]) / np.sqrt(inertia_in_lammps)
            
                        oxconf.write('%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n' % \
                          (cm[0], cm[1], cm[2], a1[0], a1[1], a1[2], a3[0], a3[1], a3[2], v[0], v[1], v[2], Lv[0], Lv[1], Lv[2]))
            
                line = lmptrj.readline()

    print("## Wrote data to '%s' / '%s'" % (configuration_file, topology_file), file=sys.stderr)
    print("## DONE", file=sys.stderr)


if __name__ == '__main__':
    main()
