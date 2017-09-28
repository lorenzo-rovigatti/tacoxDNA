import numpy as np
import math as mt




#si da in input un ndarray con le coordinate della curva, che puo essere asse del dna o uno dei due strand
def get_writhe(coordxyz):
 	numrows = len(coordxyz)
	dist=np.copy(coordxyz)
	#distn=np.copy(coordxyz)

	WR=0

	#costruiamo la matrice delle distanze
	for c in range(0,numrows): 
		ind = int( c - mt.floor(c/ float(numrows) ) * numrows		)
		ind1 =int( c+1 - mt.floor((c+1)/ float(numrows) ) * numrows	)

		dist[ind,:]=coordxyz[ind1,:]-coordxyz[ind,:]
		#distn[ind,:]=dist[ind,:]/np.sqrt(np.dot(dist[ind,:],dist[ind,:]))	
		#print ind,ind1, coordxyz[ind,:],coordxyz[ind1,:],dist[ind,:]

        #1 metodo paperKlenin Langowski 2000, uso la stessa notazione
	for i in range(1,numrows):
		#if(i%100==0):
		#	print i
		for j in range(0,i):
			#i indica il segmento i,i+1 # il segmento j,j+1
			#calcoliamo i 4 segmenti necessari
			ind_i = int(	i - mt.floor(i/ float(numrows) ) * numrows		)
			ind_i1 =int( 	i+1 - mt.floor((i+1)/ float(numrows) ) * numrows	)
			ind_j = int(	j - mt.floor(j/ float(numrows) ) * numrows		)
			ind_j1 =int( 	j+1 - mt.floor((j+1)/ float(numrows) ) * numrows	)

			r12	=	dist[ind_i,:]    			#rii1
			r34	=	dist[ind_j,:]    			#rjj1
			r13	=	coordxyz[ind_j,:] - coordxyz[ind_i,:]	#rij
			r23	=	coordxyz[ind_j,:] - coordxyz[ind_i1,:]	#ri1j
			r24	=	coordxyz[ind_j1,:] - coordxyz[ind_i1,:] #ri1j1
			r14	=	coordxyz[ind_j1,:] - coordxyz[ind_i,:]  #rij1

			#print ind_i,ind_i1,ind_j,ind_j1,r12,dist[ind_i,:],r34,r13,r23,r24,r14

			#do action only if 4 points are coplanar (from wolfram), otherwise solid angle is zero
			if(abs (np.dot( r13 , np.cross(r12,r34) )) > mt.exp(-10) ):
				n1	=	np.cross( r13 , r14 )
				n1 /= np.sqrt(np.dot( n1, n1 ))
				n2	=	np.cross( r14 , r24 )
				n2 /= np.sqrt(np.dot( n2, n2 ))
				n3	=	np.cross( r24 , r23 )
				n3 /= np.sqrt(np.dot( n3, n3 ))
				n4	=	np.cross( r23 , r13 )
				n4 /= np.sqrt(np.dot( n4, n4 ))

			
				#print ind_i,ind_j,  np.dot( r13 , np.cross(r12,r34) ),  n1,n2,n3,n4

				wr_loc = np.arcsin( np.dot(n1,n2) ) + np.arcsin( np.dot(n2,n3) ) + np.arcsin( np.dot(n3,n4) ) + np.arcsin( np.dot(n4,n1) )

				wr_loc = wr_loc * np.sign( np.dot( np.cross( r34 , r12 )  ,r13     )        )   / 4 / np.pi


				WR+= 2 * wr_loc
				#print wr_loc,WR



	return WR













			
			

