#prende in input un file x y z q0 q1 q2 q3

#in lammps il nucleotide come corpo rigido puo essere rappresentato da un quaternione. Il quaternione rappresenta una rotazione di un angolo theta attorno ad un asse 
#la rotazione identifica una matrice 3x3 applicabile a un punto per fare questa rotazione. Qui in basso ci sono le formule perconvertirle usate qui sotto
#trasformazione matrice -> quaternioni e inversa

#http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
#http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm

#La terna di vettori si ricava leggendo una per una le righe (o le colonne, mi confonde), e quelle rappresentato la terna (sarebbe come moltiplicare la matrice per (1 0 0), (0 1 0) etc)
#Probabilmente per ricavare il quaternione e' stata usata la matrice di rotazione inversa riferita alla terna (quindi la matrice trasposta), e ugualmente va considerata la stessa cosa per ritornare indietro
# ho verificato che la trasformazione all indietro coindide con quella in avanti

#nella terna del nucleotide
#a1 corrisponde al vettore congiungente backbone to base (per interazioni base-base)
#a3 corisponde al terzo vettore della terna e alla congiungente del backbone

import numpy as np
import sys, os

def exyz_to_quat (mya1, mya3):

    mya2 = np.cross(mya3, mya1)
    myquat = [1,0,0,0]

    q0sq = 0.25 * (mya1[0] + mya2[1] + mya3[2] + 1.0)
    q1sq = q0sq - 0.5 * (mya2[1] + mya3[2])
    q2sq = q0sq - 0.5 * (mya1[0] + mya3[2])
    q3sq = q0sq - 0.5 * (mya1[0] + mya2[1])

    # some component must be greater than 1/4 since they sum to 1
    # compute other components from it

    if q0sq >= 0.25:
	myquat[0] = np.sqrt(q0sq)
	myquat[1] = (mya2[2] - mya3[1]) / (4.0*myquat[0])
	myquat[2] = (mya3[0] - mya1[2]) / (4.0*myquat[0])
	myquat[3] = (mya1[1] - mya2[0]) / (4.0*myquat[0])
    elif q1sq >= 0.25:
	myquat[1] = np.sqrt(q1sq)
	myquat[0] = (mya2[2] - mya3[1]) / (4.0*myquat[1])
	myquat[2] = (mya2[0] + mya1[1]) / (4.0*myquat[1])
	myquat[3] = (mya1[2] + mya3[0]) / (4.0*myquat[1])
    elif q2sq >= 0.25:
	myquat[2] = np.sqrt(q2sq)
	myquat[0] = (mya3[0] - mya1[2]) / (4.0*myquat[2])
	myquat[1] = (mya2[0] + mya1[1]) / (4.0*myquat[2])
	myquat[3] = (mya3[1] + mya2[2]) / (4.0*myquat[2])
    elif q3sq >= 0.25:
	myquat[3] = np.sqrt(q3sq)
	myquat[0] = (mya1[1] - mya2[0]) / (4.0*myquat[3])
	myquat[1] = (mya3[0] + mya1[2]) / (4.0*myquat[3])
	myquat[2] = (mya3[1] + mya2[2]) / (4.0*myquat[3])

    norm = 1.0/np.sqrt(myquat[0]*myquat[0] + myquat[1]*myquat[1] + \
			  myquat[2]*myquat[2] + myquat[3]*myquat[3])
    myquat[0] *= norm
    myquat[1] *= norm
    myquat[2] *= norm
    myquat[3] *= norm

    return np.array([myquat[0],myquat[1],myquat[2],myquat[3]])

def quat_to_exyz (myquat,mya1, mya3):
		
	sqw = myquat[0]*myquat[0];
    	sqx = myquat[1]*myquat[1];
    	sqy = myquat[2]*myquat[2];
    	sqz = myquat[3]*myquat[3];

	invs = 1 / (sqx + sqy + sqz + sqw)
    	m00 = ( sqx - sqy - sqz + sqw)*invs ;
    	m11 = (-sqx + sqy - sqz + sqw)*invs ;
    	m22 = (-sqx - sqy + sqz + sqw)*invs ;
    
        tmp1 = myquat[1]*myquat[2];
    	tmp2 = myquat[3]*myquat[0];
    	m10 = 2.0 * (tmp1 + tmp2)*invs ;
    	m01 = 2.0 * (tmp1 - tmp2)*invs ;
    
    	tmp1 = myquat[1]*myquat[3];
    	tmp2 = myquat[2]*myquat[0];
   	m20 = 2.0 * (tmp1 - tmp2)*invs ;
    	m02 = 2.0 * (tmp1 + tmp2)*invs ;
    	tmp1 = myquat[2]*myquat[3];
    	tmp2 = myquat[1]*myquat[0];
    	m21 = 2.0 * (tmp1 + tmp2)*invs ;
    	m12 = 2.0 * (tmp1 - tmp2)*invs ; 

	mya1[0]=m00
	mya1[1]=m10
	mya1[2]=m20

	mya3[0]=m02
	mya3[1]=m12
	mya3[2]=m22

	return
	

coordxyz=np.loadtxt(sys.argv[1],float)
numrows = len(coordxyz)

nbasi=numrows/2

#non flippiamo la seconda base perche deve essere stampata in questo ordine
#ssdna1 = coordxyz[0:nbasi,0:3]
#ssdna2 = coordxyz[nbasi:nbasi*2,0:3]
#ssdna2 = np.flipud (coordxyz[nbasi:nbasi*2,:])
#a1 corrisponde a vperp nel mio programma e sarebbe il primo vettore o mya1 e darebbe la direzione verso la base dal cm
#a1_ssdna1=coordxyz[0:nbasi,3:6]
#a1_ssdna2=coordxyz[nbasi:nbasi*2,3:6]
#a3 corrisponde a vdir e mya3
#a3_ssdna1=coordxyz[0:nbasi,6:9]
#a3_ssdna2=coordxyz[nbasi:nbasi*2,6:9]
#velocita varie


min_boxx=min(coordxyz[:,0]) 
max_boxx=max(coordxyz[:,0]) 
min_boxy=min(coordxyz[:,1]) 
max_boxy=max(coordxyz[:,1]) 
min_boxz=min(coordxyz[:,2]) 
max_boxz=max(coordxyz[:,2])
max_box=max(max_boxx,max_boxy,max_boxz)
min_box=max(min_boxx,min_boxy,min_boxz)
#informazioni su tipo basi (create a caso)
basetype=np.zeros(numrows)
a1=np.zeros(3)
a3=np.zeros(3)

#for i in xrange (nbasi):
#                val_base=np.random.randint(0, 4)
#                basetype[i]=val_base
#                basetype[numrows-i-1]=(3-val_base)





out = open ("oxdna.input", "w")

out.write('t = 0\n')
out.write('b = %f %f %f\n' % (max_boxx-min_boxx+10,max_boxy-min_boxy+10,max_boxz-min_boxz+10))
out.write('E = 0. 0. 0.\n')

for i in xrange(numrows):
    quaternions=coordxyz[i,3:7]
    quat_to_exyz(quaternions,a1,a3)
    out.write('%22.15le %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le 0 0 0 0 0 0\n' \
          % (coordxyz[i,0], coordxyz[i,1], coordxyz[i,2], \
             a1[0],a1[1],a1[2],a3[0],a3[1],a3[2]  ))

    #quaternions1=exyz_to_quat(a1,a3)
    #print >> sys.stdout,quaternions,quaternions1
  

print >> sys.stdout, "## Wrote data to 'oxdna.input'"
print >> sys.stdout, "## DONE"







