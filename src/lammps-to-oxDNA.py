# prende in input un file x y z q0 q1 q2 q3

# in lammps il nucleotide come corpo rigido puo essere rappresentato da un quaternione. Il quaternione rappresenta una rotazione di un angolo theta attorno ad un asse 
# la rotazione identifica una matrice 3x3 applicabile a un punto per fare questa rotazione. Qui in basso ci sono le formule perconvertirle usate qui sotto
# trasformazione matrice -> quaternioni e inversa

# http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
# http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm

# La terna di vettori si ricava leggendo una per una le righe (o le colonne, mi confonde), e quelle rappresentato la terna (sarebbe come moltiplicare la matrice per (1 0 0), (0 1 0) etc)
# Probabilmente per ricavare il quaternione e' stata usata la matrice di rotazione inversa riferita alla terna (quindi la matrice trasposta), e ugualmente va considerata la stessa cosa per ritornare indietro
# ho verificato che la trasformazione all indietro coindide con quella in avanti

# nella terna del nucleotide
# a1 corrisponde al vettore congiungente backbone to base (per interazioni base-base)
# a3 corisponde al terzo vettore della terna e alla congiungente del backbone

import numpy as np
import sys, os
from libs import base


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
    if len(sys.argv) < 3:
        print >> sys.stderr, "Usage is %s lammps_input configuration" % sys.argv[0]
        sys.exit(1)

    with open(sys.argv[2]) as conf:
        conf.readline()
        conf.readline()
        conf.readline()
        N = int(conf.readline())
        conf.readline()
        box_min = np.array([0., 0., 0.])
        box_max = np.array([0., 0., 0.])
        for i in range(3):
            box_min[i], box_max[i] = [float(x) for x in conf.readline().split()]
        conf.readline()

        coordxyz = np.loadtxt(conf, float)
        if len(coordxyz) != N:
            print >> sys.stderr, "The number of particles specified in the headers of the configuration (%d) is different from the particle lines found therein (%d)" % (N, len(coordxyz))
            sys.exit(1)

        box = box_max - box_min
        system = base.System(box)

        current_strand = base.Strand()

        for nucleotide in coordxyz:
            cm = nucleotide[2:5]
            quaternions = nucleotide[5:9]
            a1, a3 = quat_to_exyz(quaternions)
            b = int(nucleotide[1]) - 1

            current_strand.add_nucleotide(base.Nucleotide(cm, a1, a3, b, b))

        system.add_strand(current_strand)
        basename = os.path.basename(sys.argv[2])
        system.print_lorenzo_output(basename + ".oxdna", basename + ".top")

