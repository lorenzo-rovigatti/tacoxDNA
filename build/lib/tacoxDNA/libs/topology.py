import numpy as np
import math as mt
import numpy.linalg as la

#angle between vectors (-pi,pi) with a reference direction vplane 
def py_ang(v1, v2, vplane):
	v1n = v1 / la.norm(v1)
	v2n = v2 / la.norm(v2)

	return np.arctan2 (la.norm(np.cross(v1n, v2n)) , np.dot (v1n, v2n)) * np.sign(np.dot(np.cross(v1n, v2n), vplane))


def get_twist(axis, ssdna1):
	numrows = len(axis)
	distn = np.copy(axis)
	dist = np.copy(axis)
	p = np.copy(axis)
	a = np.copy(axis)

	TW = 0

        #axis vectors
	for c in range(0, numrows):
		ind = c
		ind1 = (c + 1) % numrows

		dist[ind, :] = axis[ind1, :] - axis[ind, :]
		distn[ind, :] = dist[ind, :] / np.sqrt(np.dot(dist[ind, :], dist[ind, :]))
		# print ind,ind1, dist[ind,:],axis[ind,:],axis[ind1,:]

	# vector perpendicular to two consecutive axis vectors
	for c in range(0, numrows):
		ind_1 = (c - 1 + numrows) % numrows
		ind = c

		p[ind, :] = np.cross(dist[ind_1, :] , dist[ind, :])
		p[ind, :] /= np.sqrt(np.dot(p[ind, :] , p[ind, :]))
		# print p[ind,0], p[ind,1], p[ind,2]

	# axis to base vectors (perpendicular to axis)
	weight = 0.5  # 1 #0.5
	for c in range(0, numrows):
		a[c, :] = weight * ssdna1[c, :] + (1 - weight) * ssdna1[(c - 1 + numrows) % numrows, :] - axis[c, :]

	for c in range(0, numrows):
		proj = np.dot(a[c, :], distn[c, :])
		a[c, :] = a[c, :] - proj * distn[c, :]

	# twist angles
	for c in range(0, numrows):
		ind_1 = (c - 1 + numrows) % numrows
		ind = c

		# the angle should be computed selecting dist as the axis, so we need to choose the right order

		alpha = py_ang(a[ind_1, :] , p[ind, :] , dist[ind_1, :])
		gamma = py_ang(p[ind, :]  , a[ind, :]  , dist[ind, :])

		angle = (alpha + gamma + 4 * np.pi) % (2 * np.pi)
		# Now we have the angle in 0 - 2pi . if it exceeds pi, let's take angle-2*pi instead
		if angle > np.pi:
			angle = angle - 2 * np.pi

		# print  angle

		TW += angle / (2 * np.pi)

	return TW


#curve xyz coordinate as input
def get_writhe(coordxyz):
	numrows = len(coordxyz)
	dist = np.copy(coordxyz)
	# distn=np.copy(coordxyz)

	WR = 0

	for c in range(0, numrows): 
		ind = int(c - mt.floor(c / float(numrows)) * numrows)
		ind1 = int(c + 1 - mt.floor((c + 1) / float(numrows)) * numrows)

		dist[ind, :] = coordxyz[ind1, :] - coordxyz[ind, :]
		# distn[ind,:]=dist[ind,:]/np.sqrt(np.dot(dist[ind,:],dist[ind,:]))	
		# print ind,ind1, coordxyz[ind,:],coordxyz[ind1,:],dist[ind,:]

	for i in range(1, numrows):
		for j in range(0, i):
			ind_i = int(i - mt.floor(i / float(numrows)) * numrows)
			ind_i1 = int(i + 1 - mt.floor((i + 1) / float(numrows)) * numrows)
			ind_j = int(j - mt.floor(j / float(numrows)) * numrows)
			ind_j1 = int(j + 1 - mt.floor((j + 1) / float(numrows)) * numrows)

			r12	 = 	dist[ind_i, :]  # rii1
			r34	 = 	dist[ind_j, :]  # rjj1
			r13	 = 	coordxyz[ind_j, :] - coordxyz[ind_i, :]  # rij
			r23	 = 	coordxyz[ind_j, :] - coordxyz[ind_i1, :]  # ri1j
			r24	 = 	coordxyz[ind_j1, :] - coordxyz[ind_i1, :]  # ri1j1
			r14	 = 	coordxyz[ind_j1, :] - coordxyz[ind_i, :]  # rij1

			# print ind_i,ind_i1,ind_j,ind_j1,r12,dist[ind_i,:],r34,r13,r23,r24,r14

			# do action only if 4 points are coplanar (from wolfram), otherwise solid angle is zero
			# if(ind_i==ind_j+1):
			# 	if(np.dot( r13 , np.cross(r12,r34) )> 5*mt.exp(-17)):
			# 		print np.dot( r13 , np.cross(r12,r34) ),ind_i,ind_j
			if(abs (np.dot(r13 , np.cross(r12, r34))) > 5 * mt.exp(-17)):
				n1	 = 	np.cross(r13 , r14)
				n1 /= np.sqrt(np.dot(n1, n1))
				n2	 = 	np.cross(r14 , r24)
				n2 /= np.sqrt(np.dot(n2, n2))
				n3	 = 	np.cross(r24 , r23)
				n3 /= np.sqrt(np.dot(n3, n3))
				n4	 = 	np.cross(r23 , r13)
				n4 /= np.sqrt(np.dot(n4, n4))
			
				# print ind_i,ind_j,  np.dot( r13 , np.cross(r12,r34) ),  n1,n2,n3,n4

				wr_loc = np.arcsin(np.dot(n1, n2)) + np.arcsin(np.dot(n2, n3)) + np.arcsin(np.dot(n3, n4)) + np.arcsin(np.dot(n4, n1))

				wr_loc = wr_loc * np.sign(np.dot(np.cross(r34 , r12)  , r13)) / 4 / np.pi

				WR += 2 * wr_loc
				# print wr_loc,WR

	return WR
