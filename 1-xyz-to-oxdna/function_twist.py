import numpy as np
import math as mt
import numpy.linalg as la
#versione del 31 MARZO 2017

#angolo definito dando una normale fissata del piano
def py_ang(v1, v2,vplane):
	v1n = v1 / la.norm( v1 )
	v2n = v2 / la.norm( v2 )


	#arctan con questi argomenti da angolo tra 0 e pi, e la normale definisce la direzione. Se vuoi fissare la direzione, il segno lo decide quale delle due normali scegli
	return np.arctan2 (  la.norm(np.cross(v1n, v2n)) , np.dot (v1n, v2n) ) * np.sign( np.dot( np.cross(v1n,v2n), vplane) )





#in questo caso diamo l asse e uno strand del dna, e abbiamo tutti i dati, seguiamo metodo 3 del paper di langowski
def get_twist(axis,ssdna1):
	numrows = len(axis)
	distn=np.copy(axis)
	dist=np.copy(axis)
	p=np.copy(axis)
	a=np.copy(axis)

	TW=0

	#calcoliamo vettori distanza tra elementi successivi asse
	for c in range(0,numrows):
		ind = c
		ind1 = (c+1)%numrows

		dist[ind,:]=axis[ind1,:]-axis[ind,:]
		distn[ind,:]=dist[ind,:]/np.sqrt(np.dot(dist[ind,:],dist[ind,:]))
		#print ind,ind1, dist[ind,:],axis[ind,:],axis[ind1,:]

	#calcoliamo vettori perpendicolare a due segmenti successivi
	for c in range(0,numrows):
		ind_1 =  (c-1+numrows)%numrows
		ind = c

		p[ind,:] = np.cross( dist[ind_1,:] , dist[ind,:] )
		p[ind,:]/= np.sqrt(np.dot( p[ind,:] , p[ind,:] ))
		#print p[ind,0], p[ind,1], p[ind,2]


	#calcoliamo vettori asse-nucleotide e ortogonalizziamo
	weight=0.5 #1 #0.5
	for c in range(0,numrows):
		a[c,:]=weight*ssdna1[c,:] + (1-weight)*ssdna1[(c-1+numrows)%numrows,:]-axis[c,:]

	for c in range(0,numrows):
		proj=np.dot(a[c,:],distn[c,:])
		a[c,:]=a[c,:]-proj*distn[c,:]

	#calcoliamo gli angoli di twist separatamente
	for c in range(0,numrows):
		ind_1 = (c-1+numrows)%numrows
		ind = c

		#the angle should be computed selecting dist as the axis, so we need to choose the right order

		alpha = py_ang( a[ind_1,:] , p[ind,:] , dist[ind_1,:])
		gamma = py_ang(	p[ind,:]  , a[ind,:]  , dist[ind,:])

		#vogliamo che resti tra -pi e pi

		angle= (alpha+gamma + 4*np.pi) % (2*np.pi)
		# Now we have the angle in 0 - 2pi . if it exceeds pi, let's take angle-2*pi instead
		if angle >np.pi:
			angle= angle - 2*np.pi

		#print  angle

		TW += angle / (2*np.pi)


	return TW


















