#!/usr/bin/env python3

import numpy as np
import sys, os
from libs import base
from libs import reader_lammps_init 
from libs.constants import mass_in_lammps, inertia_in_lammps, number_oxdna_to_lammps

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


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage is %s lammps_init_file" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    conf = reader_lammps_init.Lammps_parser(sys.argv[1])
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
            b = number_oxdna_to_lammps[conf.bases[i]-1] 

            v = np.array(conf.v[i,:]) * np.sqrt(mass_in_lammps)
            Lv = np.array(conf.Lv[i,:]) / np.sqrt(inertia_in_lammps)

            strands[conf.strand[i]-1].add_nucleotide(base.Nucleotide(cm, a1, a3, b, b, v, Lv))

            # close strand 
            next_bond=conf.bonds[i][1]
            if next_bond!=-1 and next_bond!=i+1:
                if conf.strand[i]!=conf.strand[next_bond]:
                    print("Wrong bond arising between two different strands", file=sys.stderr)
                else:
                    strands[conf.strand[i]-1].make_circular()


    for i in range(conf.nstrands):
        system.add_strand(strands[i])

    basename = os.path.basename(sys.argv[1])
    topology_file = basename + ".top"
    configuration_file = basename + ".oxdna"
    system.print_lorenzo_output(configuration_file, topology_file)

    if len(sys.argv) == 3:

        f = open(sys.argv[2],'r')
        t = []

        line = f.readline()

        while line:

            if line.startswith('ITEM: TIMESTEP'):
                t = int(f.readline())
                line = f.readline()

            if line.startswith('ITEM: NUMBER OF ATOMS') and t==0:
                conf.natoms = int(f.readline())
                line = f.readline()

            if line.startswith('ITEM: BOX BOUNDS') and t==0:
                line = f.readline()
                xlo, xhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                conf.Lx = xhi - xlo
                line = f.readline()
                ylo, yhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                conf.Ly = yhi - ylo
                line = f.readline()
                zlo, zhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                conf.Lz = zhi - zlo

            if line.startswith('ITEM: ATOMS'):
                line = f.readline()
                N = conf.natoms

                # converting LAMMPS data into native oxDNA data format
                for n in range(N):

                    index = int(line.split()[0])-1
                    # Position
                    cm = np.float32(line.split()[3:6])
                    conf.xyz[index,:] = cm
                    # Velocity 
                    velocity = np.float32(line.split()[13:16])
                    conf.v[index,:] = velocity
                    # Quaternions 
                    quaternions = np.float32(line.split()[9:13])
                    conf.ellipsoids[index,:] = quaternions
                    # Angular momentum 
                    angmom = np.float32(line.split()[16:19]) 
                    conf.Lv[index,:] = angmom

                    line = f.readline()

                # write oxDNA data to file

                cf = open(configuration_file,'a')
                # header goes here
                for n in range(N):

                    cm = conf.xyz[n,:]
                    quaternions = conf.ellipsoids[n,:]
                    a1, a3 = quat_to_exyz(quaternions)
                    v = np.array(conf.v[n,:]) * np.sqrt(mass_in_lammps)
                    Lv = np.array(conf.Lv[n,:]) / np.sqrt(inertia_in_lammps)

                    print(cm, quaternions, a1,a3,v,Lv)

                    cf.write('%le %le %le\n' % (cm[0],cm[1],cm[2]))

            line = f.readline()

        f.close()
        cf.close()

    print("## Wrote data to '%s' / '%s'" % (configuration_file, topology_file), file=sys.stderr)
    print("## DONE", file=sys.stderr)
