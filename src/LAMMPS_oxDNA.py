import numpy as np
import sys, os
from libs import base
from libs import reader_lammps_init 

number_oxdna_to_lammps = {0 : 0, 1 : 2, 2 : 1, 3 : 3}


#rules to convert from a1 a3 vectors to quaternions based on
# http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
# http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm

def exyz_to_quat (mya1, mya3):

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


def quat_to_exyz (myquat):
    sqw = myquat[0] * myquat[0];
    sqx = myquat[1] * myquat[1];
    sqy = myquat[2] * myquat[2];
    sqz = myquat[3] * myquat[3];

    invs = 1 / (sqx + sqy + sqz + sqw)
    m00 = (sqx - sqy - sqz + sqw) * invs ;
    m11 = (-sqx + sqy - sqz + sqw) * invs ;
    m22 = (-sqx - sqy + sqz + sqw) * invs ;
    
    tmp1 = myquat[1] * myquat[2];
    tmp2 = myquat[3] * myquat[0];
    m10 = 2.0 * (tmp1 + tmp2) * invs ;
    m01 = 2.0 * (tmp1 - tmp2) * invs ;

    tmp1 = myquat[1] * myquat[3];
    tmp2 = myquat[2] * myquat[0];
    m20 = 2.0 * (tmp1 - tmp2) * invs ;
    m02 = 2.0 * (tmp1 + tmp2) * invs ;
    tmp1 = myquat[2] * myquat[3];
    tmp2 = myquat[1] * myquat[0];
    m21 = 2.0 * (tmp1 + tmp2) * invs ;
    m12 = 2.0 * (tmp1 - tmp2) * invs ; 

    mya1 = np.array([m00, m10, m20])
    mya3 = np.array([m02, m12, m22])

    return mya1, mya3


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print >> sys.stderr, "Usage is %s lammps_init_file" % sys.argv[0]
        sys.exit(1)

    conf=reader_lammps_init.Lammps_parser(sys.argv[1])
    N = conf.natoms 
    box = np.array([0, 0., 0.])
    box[0]=conf.Lx
    box[1]=conf.Ly
    box[2]=conf.Lz

    system = base.System(box)


    strands=[]
    for i in range(conf.nstrands):
            strands.append(base.Strand())


    for i in range(N):
            cm = conf.xyz[i,:]
            quaternions = conf.ellipsoids[i,:]
            a1, a3 = quat_to_exyz(quaternions)
            b = number_oxdna_to_lammps[conf.bases[i]-1] 

            v = conf.v[i,:]
            Lv = conf.Lv[i,:]

            strands[conf.strand[i]-1].add_nucleotide(base.Nucleotide(cm, a1, a3, b, b,v,Lv))

            #close strand 
            next_bond=conf.bonds[i][1]
            if next_bond!=-1 and next_bond!=i+1:
                if conf.strand[i]!=conf.strand[next_bond]:
                    print >> sys.stderr, "Wrong bond arising between two different strands"
                else:
                    strands[conf.strand[i]-1].make_circular()


    for i in range(conf.nstrands):
        system.add_strand(strands[i])

    basename = os.path.basename(sys.argv[1])
    system.print_lorenzo_output(basename + ".oxdna", basename + ".top")

