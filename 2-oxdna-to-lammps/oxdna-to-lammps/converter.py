#convertiamo oxdna to lammps starter
#per ora facciamo solo il caso con 2 strand uguali
#testa mettendo il bo cubico invece che rettangolare

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
max_boxx=min_boxx+800   #   max(coordxyz[:,0]) 
min_boxy=min(coordxyz[:,1]) 
max_boxy=min_boxy+800 #  max(coordxyz[:,1]) 
min_boxz=min(coordxyz[:,2]) 
max_boxz=min_boxz+800 #    max(coordxyz[:,2])
max_box=max(max_boxx,max_boxy,max_boxz)
min_box=max(min_boxx,min_boxy,min_boxz)
#informazioni su tipo basi (create a caso)
basetype=np.zeros(numrows)


for i in xrange (nbasi):
                val_base=np.random.randint(0, 4)
                basetype[i]=val_base
                basetype[numrows-i-1]=(3-val_base)

print >> sys.stdout,numrows,np.absolute(max_boxx+5-min_boxx-5),np.absolute(max_boxy+5-min_boxy-5),np.absolute(max_boxz+5-min_boxz-5)
print >> sys.stdout, (4./3.)*numrows*np.pi/np.absolute(max_boxx+5-min_boxx-5)/np.absolute(max_boxy+5-min_boxy-5)/np.absolute(max_boxz+5-min_boxz-5)


out = open ("test_oxdna.input", "w")

out.write('# LAMMPS data file\n')
out.write('%d atoms\n' % numrows)
out.write('%d ellipsoids\n' % numrows)
out.write('%d bonds\n' % numrows)
out.write('\n')
out.write('4 atom types\n')
out.write('1 bond types\n')
out.write('\n')
out.write('# System size\n')
out.write('%f %f xlo xhi\n' % (min_boxx-5,max_boxx+5))
out.write('%f %f ylo yhi\n' % (min_boxy-5,max_boxy+5))
out.write('%f %f zlo zhi\n' % (min_boxz-5,max_boxz+5))

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

for i in xrange(numrows):
    out.write('%d %d %22.15le %22.15le %22.15le %d 1 1\n' \
          % (i+1, basetype[i]+1, \
             coordxyz[i,0], coordxyz[i,1], coordxyz[i,2], \
             int(i/nbasi)+1))

out.write('\n')
out.write('# Atom-ID, translational, rotational velocity\n')
out.write('Velocities\n')
out.write('\n')

for i in xrange(numrows):
    #out.write("%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n" \
    #      % (i+1,0.0,0.0,0.0,0.0,0.0,0.0))
    out.write("%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n" \
          % (i+1,coordxyz[i,9],coordxyz[i,10],coordxyz[i,11],coordxyz[i,12],coordxyz[i,13],coordxyz[i,14]))

out.write('\n')
out.write('# Atom-ID, shape, quaternion\n')
out.write('Ellipsoids\n')
out.write('\n')

for i in xrange(numrows):
    quaternions=exyz_to_quat(coordxyz[i,3:6],coordxyz[i,6:9])
    out.write(\
    "%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n"  \
      % (i+1,1.1739845031423408,1.1739845031423408,1.1739845031423408, \
    quaternions[0],quaternions[1], quaternions[2],quaternions[3]))

out.write('\n')
out.write('# Bond topology\n')
out.write('Bonds\n')
out.write('\n')

for i in xrange(numrows):
    if i+1==nbasi or i+1==2*nbasi:
        out.write("%d  %d  %d  %d\n" % (i+1,1,i+1,i+2-nbasi))
    else:
        out.write("%d  %d  %d  %d\n" % (i+1,1,i+1,i+2))

out.close()

print >> sys.stdout, "## Wrote data to 'data.oxdna'"
print >> sys.stdout, "## DONE"







